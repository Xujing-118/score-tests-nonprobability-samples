rm(list = ls())

library(dplyr)
library(MASS)

load("D:/paper/MAR/espacov/sampA.RData")
load("D:/paper/MAR/espacov/sampB.RData")

SA <- as.data.frame(sampA)
SB <- as.data.frame(sampB)

prop_vars <- c("age1", "age2", "age3", "educ", "gender", "assessGovernment", "cost", "vaccine", "socialStat")
outcome_vars <- c(prop_vars, "health")

SA <- SA[, c(outcome_vars, "y")]
SB <- SB[, c(outcome_vars, "weight")]

names(SA)[names(SA) == "y"] <- "Y"
names(SB)[names(SB) == "weight"] <- "weight_raw"

SA$R <- 1L
SA$weight_raw <- 1

SB$Y <- NA_real_
SB$R <- 0L

SA <- SA[, c("Y", outcome_vars, "R", "weight_raw")]
SB <- SB[, c("Y", outcome_vars, "R", "weight_raw")]

dat <- bind_rows(SA, SB)

idxA <- which(dat$R == 1)
idxB <- which(dat$R == 0)

YA <- dat$Y[idxA]

dB <- as.numeric(dat$weight_raw[idxB])
target_N <- sum(dB)

rhs <- function(vars) {
  if (length(vars) == 0) return("1")
  paste(vars, collapse = " + ")
}

make_X <- function(dat, vars, reference_index) {
  X <- model.matrix(as.formula(paste("~", rhs(vars))), data = dat)
  qx <- qr(X[reference_index, , drop = FALSE])
  keep <- sort(qx$pivot[seq_len(qx$rank)])
  X[, keep, drop = FALSE]
}

softplus <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

safe_inv <- function(M) {
  tryCatch(solve(M), error = function(e) MASS::ginv(M))
}

quad_form <- function(v, M_inv) {
  as.numeric(t(v) %*% M_inv %*% v)
}

# Participation model under H0

fit_beta <- function(ZA, ZB, dB) {
  scale_N <- sum(dB)
  
  fn <- function(beta) {
    etaA <- drop(ZA %*% beta)
    etaB <- drop(ZB %*% beta)
    -(sum(etaA) - sum(dB * softplus(etaB))) / scale_N
  }
  
  gr <- function(beta) {
    piB <- plogis(drop(ZB %*% beta))
    -(colSums(ZA) - drop(crossprod(ZB, dB * piB))) / scale_N
  }
  
  beta0 <- rep(0, ncol(ZA))
  p0 <- nrow(ZA) / sum(dB)
  beta0[1] <- qlogis(p0)
  
  fit <- optim(par = beta0, fn = fn, gr = gr, method = "BFGS", control = list(maxit = 10000, reltol = 1e-11))
  
  fit$par
}

# Probability-sample variance component

calc_HT_poisson <- function(phiB, dB, target_N) {
  piB <- 1 / dB
  sum((1 - piB) * (phiB / piB)^2) / target_N
}

# Outcome model:health is included only in the outcome model

W_all <- make_X(dat, outcome_vars, reference_index = idxA)

WA <- W_all[idxA, , drop = FALSE]
WB <- W_all[idxB, , drop = FALSE]

fit_y <- glm.fit(x = WA, y = YA, family = binomial(link = "logit"))
coef_y <- fit_y$coefficients

muA <- plogis(drop(WA %*% coef_y))
muB <- plogis(drop(WB %*% coef_y))

dmuA <- muA * (1 - muA)

score_yA <- WA * as.numeric(YA - muA)

# Score tests

run_one <- function(prop_vars, model_name) {
  
  Z_all <- make_X(dat, prop_vars, reference_index = idxA)
  
  ZA <- Z_all[idxA, , drop = FALSE]
  ZB <- Z_all[idxB, , drop = FALSE]
  
  beta_hat <- fit_beta(ZA, ZB, dB)
  
  prA <- plogis(drop(ZA %*% beta_hat))
  prB <- plogis(drop(ZB %*% beta_hat))
  
  eps <- 1e-12
  prA <- pmin(pmax(prA, eps), 1 - eps)
  prB <- pmin(pmax(prB, eps), 1 - eps)
  
  N_hat_A <- sum(1 / prA)
  
  A_hat <- crossprod(ZA, ZA * as.numeric(1 - prA)) / N_hat_A
  invA <- safe_inv(A_hat)
  
  B_hat <- crossprod(score_yA) / N_hat_A
  invB <- safe_inv(B_hat)

  S1 <- sum(YA * (1 - prA) + prA * muA) - sum(dB * prB * muB)  # S1

  A1_hat <- drop(crossprod(ZA, (1 - prA) * YA)) / N_hat_A
  B1_hat <- drop(crossprod(WA, (1 - prA) * dmuA)) / N_hat_A
  
  A2_hat <- sum((1 - prA)^2 * YA^2) / N_hat_A
  B2_hat <- sum(prA * (1 - prA) * muA^2) / N_hat_A
  
  sigma11 <- A2_hat + B2_hat - quad_form(A1_hat, invA) - quad_form(B1_hat, invB)
  
  phi1_B <- prB * (muB - drop(ZB %*% (invA %*% A1_hat)))
  
  sigma12 <- calc_HT_poisson(phi1_B, dB, target_N)
  sigma1_total <- sigma11 + sigma12
  
  Z1 <- S1 / sqrt(N_hat_A * sigma1_total)
  P1 <- 2 * pnorm(abs(Z1), lower.tail = FALSE)
  
  S2 <- sum(muA) - sum(dB * prB * muB)  # S2
  
  A3_hat <- drop(crossprod(ZA, (1 - prA) * muA)) / N_hat_A
  A4_hat <- sum((1 - prA) * muA^2) / N_hat_A
  
  sigma21 <- A4_hat - quad_form(A3_hat, invA)
  
  phi2_B <- prB * (muB - drop(ZB %*% (invA %*% A3_hat)))
  
  sigma22 <- calc_HT_poisson(phi2_B, dB, target_N)
  sigma2_total <- sigma21 + sigma22
  
  Z2 <- S2 / sqrt(N_hat_A * sigma2_total)
  P2 <- 2 * pnorm(abs(Z2), lower.tail = FALSE)
  
  data.frame(
    Model = model_name,
    N_hat_A = N_hat_A,
    S1 = S1,
    sigma11 = sigma11,
    sigma12 = sigma12,
    sigma1_total = sigma1_total,
    SE1 = sqrt(N_hat_A * sigma1_total),
    Z1 = Z1,
    P1 = P1,
    S2 = S2,
    sigma21 = sigma21,
    sigma22 = sigma22,
    sigma2_total = sigma2_total,
    SE2 = sqrt(N_hat_A * sigma2_total),
    Z2 = Z2,
    P2 = P2
  )
}


result <- run_one(prop_vars = prop_vars, model_name = "Full")
result_print <- result

num_cols <- vapply(result_print, is.numeric, logical(1))
result_print[num_cols] <- lapply(result_print[num_cols], function(x) formatC(x, format = "f", digits = 8))
print(result_print, row.names = FALSE)

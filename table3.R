rm(list = ls())

library(readxl)
library(dplyr)
library(MASS)


SA <-file.path("data", "SA1.xlsx")
SB <-file.path("data", "SB1.xlsx")

Y_name <- "SOCTRUST"
weight_name <- "WEIGHT"
cov_names <- c("AGECAT", "GENDER", "EDUC5")

# 2021 ACS
target_N <- 253554816

SA_sub <- SA[, c(Y_name, cov_names)]
SB_sub <- SB[, c(cov_names, weight_name)]

names(SA_sub)[names(SA_sub) == Y_name] <- "Y"
names(SB_sub)[names(SB_sub) == weight_name] <- "weight_raw"

SA_sub$Y <- as.numeric(SA_sub$Y)
SA_sub$R <- 1L
SA_sub$weight_raw <- 1

SB_sub$Y <- NA_real_
SB_sub$R <- 0L
SB_sub$weight_raw <- as.numeric(SB_sub$weight_raw)

SA_sub <- SA_sub[, c("Y", cov_names, "R", "weight_raw")]
SB_sub <- SB_sub[, c("Y", cov_names, "R", "weight_raw")]

dat <- bind_rows(SA_sub, SB_sub)

for (v in cov_names) {
  dat[[v]] <- factor(as.character(dat[[v]]))
}


idxA <- which(dat$R == 1)
idxB <- which(dat$R == 0)

YA <- dat$Y[idxA]

dB <- dat$weight_raw[idxB] *target_N / sum(dat$weight_raw[idxB])

cat("sum SB weight：", sum(dB), "\n")


models <- list(
  #M1 = character(0),
  #M2 = "AGECAT",
  #M3 = "GENDER",
  #M4 = "EDUC5",
  #M5 = c("AGECAT", "GENDER"),
  M6 = c("AGECAT", "EDUC5"),
  #M7 = c("GENDER", "EDUC5"),
  M8 = c("AGECAT", "GENDER", "EDUC5")
)


rhs <- function(vars) {
  if (length(vars) == 0) return("1")
  paste(vars, collapse = " + ")
}

make_X <- function(dat, vars, reference_index) {
  
  X <- model.matrix(
    as.formula(paste("~", rhs(vars))),
    data = dat
  )
  
  qx <- qr(X[reference_index, , drop = FALSE])
  
  keep <- sort( qx$pivot[seq_len(qx$rank)] )
  
  X[, keep, drop = FALSE]
}

softplus <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

safe_inv <- function(M) {
  tryCatch(
    solve(M),
    error = function(e) MASS::ginv(M)
  )
}

quad_form <- function(v, M_inv) {
  as.numeric(t(v) %*% M_inv %*% v)
}

#  beta

fit_beta <- function(ZA, ZB, dB) {
  
  scale_N <- sum(dB)
  #scale_N <- nrow(ZA)
  fn <- function(beta) {
    etaA <- drop(ZA %*% beta)
    etaB <- drop(ZB %*% beta)
    -(sum(etaA) -sum(dB * softplus(etaB))) / scale_N
  }
  
  gr <- function(beta) {
    piB <- plogis(drop(ZB %*% beta))
    -(colSums(ZA) - drop(crossprod(ZB, dB * piB))) / scale_N
  }
  
  beta0 <- rep(0, ncol(ZA))
  
  p0 <- nrow(ZA) / sum(dB)
  p0 <- min(max(p0, 1e-12), 1 - 1e-12)
  
  beta0[1] <- qlogis(p0)
  
  fit <- optim(
    par = beta0,
    fn = fn,
    gr = gr,
    method = "BFGS",
    control = list(
      maxit = 10000,
      reltol = 1e-11
    )
  )
  
  if (fit$convergence != 0) {
    warning("Participation model not fully converged.")
  }
  
  fit$par
}


# Poisson working design variance

calc_HT_poisson <- function(phiB, dB, target_N) {
  
  piB <- 1 / dB
  if (any(piB > 1)) {
    stop("probability over 1")
  }
  sum((1 - piB) *(phiB / piB)^2 ) / target_N
}


# Outcome model

W_all <- make_X(dat,cov_names,reference_index = idxA)

WA <- W_all[idxA, , drop = FALSE]
WB <- W_all[idxB, , drop = FALSE]

fit_y <- glm.fit(x = WA,y = YA,family = binomial(link = "logit"))

if (!fit_y$converged ||
    any(!is.finite(fit_y$coefficients))) {
  stop("Outcome model fit fail")
}

coef_y <- fit_y$coefficients

muA <- plogis(drop(WA %*% coef_y))
muB <- plogis(drop(WB %*% coef_y))

dmuA <- muA * (1 - muA)

# Bernoulli logistic outcome score
score_yA <- WA * as.numeric(YA - muA)

run_one <- function(prop_vars, model_name) {
  
  Z_all <- make_X(dat,prop_vars,reference_index = idxA)
  
  ZA <- Z_all[idxA, , drop = FALSE]
  ZB <- Z_all[idxB, , drop = FALSE]
  
  beta_hat <- fit_beta(ZA, ZB, dB)
  
  prA <- plogis(drop(ZA %*% beta_hat))
  prB <- plogis(drop(ZB %*% beta_hat))
  
  eps <- 1e-12
  
  prA <- pmin(
    pmax(prA, eps),
    1 - eps
  )
  
  prB <- pmin(
    pmax(prB, eps),
    1 - eps
  )
  
  
  N_hat_A <- sum(1 / prA)
  
  A_hat <- crossprod(ZA,ZA * as.numeric(1 - prA)) / N_hat_A
  invA <- safe_inv(A_hat)
  B_hat <- crossprod(score_yA) / N_hat_A
  invB <- safe_inv(B_hat)
  
  
  # S1
  S1 <- sum(YA * (1 - prA) +prA * muA) - sum(dB * prB * muB)
  
  A1_hat <- drop(crossprod(ZA,(1 - prA) * YA)) / N_hat_A

  # Bernoulli-logistic outcome model:
  # d mu(x; xi) / d xi = mu(1-mu) W
  B1_hat <- drop( crossprod(WA,(1 - prA) * dmuA) ) / N_hat_A
  A2_hat <- sum((1 - prA)^2 * YA^2 ) / N_hat_A
  B2_hat <- sum( prA * (1 - prA) * muA^2 ) / N_hat_A
  sigma11 <-A2_hat +B2_hat -quad_form(A1_hat, invA) -quad_form(B1_hat, invB)
  
  #sigma11  become slightly negative, we will apply below
  #proj_A1 <- drop(ZA %*% (invA %*% A1_hat) )
  #proj_B1 <- drop(score_yA %*% (invB %*% B1_hat))
  #u_R1 <- (1 - prA) * (YA - proj_A1) - proj_B1
  #sigma11 <- sum(u_R1^2 + prA * (1 - prA) * (proj_A1 - muA)^2 ) / N_hat_A
  
  phi1_B <- prB * (muB -drop(ZB %*% (invA %*% A1_hat)))
  
  sigma12 <- calc_HT_poisson( phiB = phi1_B,dB = dB,target_N = sum(dB) )
  
  sigma1_total <- sigma11 + sigma12
  
  if (is.finite(sigma1_total) && sigma1_total > 0) {
    Z1 <- S1 / sqrt(N_hat_A * sigma1_total)
    P1 <- 2 * pnorm(abs(Z1), lower.tail = FALSE)
  } else {
    warning(
      model_name,
      ": total analytic variance for Z1 is nonpositive."
    )
    Z1 <- NA_real_
    P1 <- NA_real_
  }
  
  # S2
  S2 <- sum(muA) -sum(dB * prB * muB)
  
  
  A3_hat <- drop( crossprod(ZA,(1 - prA) * muA)) / N_hat_A
  A4_hat <- sum((1 - prA) * muA^2) / N_hat_A
  
  sigma21 <- A4_hat -quad_form(A3_hat, invA)
  
  phi2_B <- prB * ( muB - drop(ZB %*% (invA %*% A3_hat)))
  
  sigma22 <- calc_HT_poisson(   phiB = phi2_B, dB = dB, target_N = sum(dB)  )
  
  sigma2_total <- sigma21 + sigma22
  
  if (is.finite(sigma2_total) && sigma2_total > 0) {
    Z2 <- S2 / sqrt(N_hat_A * sigma2_total)
    P2 <- 2 * pnorm(abs(Z2), lower.tail = FALSE)
  } else {
    warning(
      model_name,
      ": total analytic variance for Z2 is nonpositive."
    )
    Z2 <- NA_real_
    P2 <- NA_real_
  }
  
  
  data.frame(
    Model = model_name,
    AGECAT = as.integer("AGECAT" %in% prop_vars),
    GENDER = as.integer("GENDER" %in% prop_vars),
    EDUC5 = as.integer("EDUC5" %in% prop_vars),
    N_hat_A = N_hat_A,
    Z1 = Z1,
    P1 = P1,
    Z2 = Z2,
    P2 = P2
  )
}

result <- bind_rows(
  lapply(names(models), function(m) {
    cat("Running", m, "...\n")
    run_one(models[[m]], m)
  })
)

result_print <- result

num_cols <- vapply(result_print,is.numeric,logical(1))

result_print[num_cols] <- lapply(
  result_print[num_cols],
  function(x) formatC(x, format = "f", digits = 10)
)

print(result_print, row.names = FALSE)
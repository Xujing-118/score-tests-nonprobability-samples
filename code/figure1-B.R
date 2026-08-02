rm(list = ls())
library(nleqslv)
library(Matrix)

set.seed(2026)

N <- 100000
nA <- 500 #500,1000,2000
nB <- 500 #500,1000,2000
sim_times <- 2000

# ============================================================
# Case B: sampled outcome distribution
# Y | X, R = 1 ~ N(-X + X^2, exp(2))
# ============================================================

x <- rnorm(N, mean = 1, sd = 1)

X_pop <- cbind(1, x, x^2)   # outcome model: (1, x, x^2)
Z_pop <- cbind(1, x)        # inclusion model: (1, x)

mu_A <- -1 * x + 1 * x^2
sigma2_true <- exp(2)
sigma_true <- sqrt(sigma2_true)

# ============================================================
# Outcome model MLE
# Y | X, R = 1 ~ N(xi1*x + xi2*x^2, exp(xi3))
# ============================================================

xiloglik <- function(xi, x, y, d) {
  res <- y - xi[1] * x - xi[2] * x^2
  sum(d * (xi[3] + res^2 * exp(-xi[3])))
}

max.xiloglik <- function(xi, x, y, d) {
  fit <- nlminb(xi, xiloglik, x = x, y = y, d = d)
  fit$par
}

# ============================================================
# H0: gamma = 0
# ============================================================

gamma_val <- 0
set.seed(2026)
# c_A(x; gamma, zeta) = log E[exp(-gamma Y) | X, R=1]
# For normal sampled outcome distribution:
# c_A = -gamma * mu_A + 0.5 * gamma^2 * sigma2
c_gamma <- -gamma_val * mu_A + 0.5 * gamma_val^2 * sigma2_true

beta0_opt <- uniroot(function(be0) {
  sum(plogis(be0 + 0.5 * x - c_gamma)) - nA
}, c(-20, 20))$root

pi_x <- plogis(beta0_opt + 0.5 * x - c_gamma)

# Generate finite population under Case B
# When gamma = 0, the two components have the same mean mu_A,
# but we keep this form to match the Case B DGP.
G <- rbinom(N, 1, pi_x)

Y <- rnorm(
  N,
  mean = ifelse(G == 1, mu_A, mu_A - gamma_val * sigma2_true),
  sd = sigma_true
)

pi_A_pop <- plogis(beta0_opt + 0.5 * x + gamma_val * Y)

cat(sprintf("H0: gamma = %.2f\n", gamma_val))
cat(sprintf(
  "beta0 = %.4f, sum(pi_x) = %.2f, sum(pi_A_pop) = %.2f\n",
  beta0_opt, sum(pi_x), sum(pi_A_pop)
))

# ============================================================
# Reference probability sample S_B
# ============================================================

temp_val <- x + 0.03 * Y
c_val <- (max(temp_val) - 50 * min(temp_val)) / 49

if (c_val + min(temp_val) <= 0) {
  c_val <- abs(min(temp_val)) + 0.1
}

z_size <- c_val + temp_val
pi_B_pop <- nB * z_size / sum(z_size)
pi_B_pop <- pmin(pi_B_pop, 0.99)

cat(sprintf("sum(pi_B_pop) = %.2f\n", sum(pi_B_pop)))
cat(sprintf("pi_B capped count = %d\n", sum(pi_B_pop >= 0.99)))

# ============================================================
# Monte Carlo simulation under H0
# ============================================================

p_values_s1 <- rep(NA_real_, sim_times)
p_values_s2 <- rep(NA_real_, sim_times)

count_neg_s1 <- 0
count_neg_s2 <- 0
count_fail_beta <- 0

cat(sprintf("Start %d simulations...\n", sim_times))

for (b in 1:sim_times) {
  
  if (b %% 500 == 0) {
    cat(sprintf("  %d / %d ...\n", b, sim_times))
  }
  
  # ------------------------------------------------------------
  # Draw S_A and S_B
  # ------------------------------------------------------------
  
  R_A <- rbinom(N, 1, pi_A_pop)
  R_B <- rbinom(N, 1, pi_B_pop)
  
  SA_x <- X_pop[R_A == 1, , drop = FALSE]
  SA_y <- Y[R_A == 1]
  
  SB_x <- X_pop[R_B == 1, , drop = FALSE]
  SB_pi <- pi_B_pop[R_B == 1]
  SB_d <- 1 / SB_pi
  
  if (length(SA_y) < 10 || nrow(SB_x) < 10) {
    next
  }
  
  # ------------------------------------------------------------
  # Estimate sampled outcome model using S_A
  # ------------------------------------------------------------
  
  xi_init <- c(-1, 1, 2)
  
  xi_hat <- max.xiloglik(
    xi_init,
    SA_x[, 2],
    SA_y,
    rep(1, length(SA_y))
  )
  
  SA_mu <- as.numeric(xi_hat[1] * SA_x[, 2] + xi_hat[2] * SA_x[, 3])
  SB_mu <- as.numeric(xi_hat[1] * SB_x[, 2] + xi_hat[2] * SB_x[, 3])
  
  SA_sigma2 <- exp(xi_hat[3])
  
  # ------------------------------------------------------------
  # Estimate beta under H0
  # ------------------------------------------------------------
  
  SA_z <- SA_x[, 1:2, drop = FALSE]
  SB_z <- SB_x[, 1:2, drop = FALSE]
  
  score_eqn <- function(beta) {
    pi_B_est <- plogis(as.numeric(SB_z %*% beta))
    
    sum_SA <- colSums(SA_z)
    sum_SB <- colSums(SB_z * as.numeric(SB_d * pi_B_est))
    
    sum_SA - sum_SB
  }
  
  sol <- nleqslv(rep(0, 2), score_eqn)
  
  if (sol$termcd > 2) {
    count_fail_beta <- count_fail_beta + 1
    next
  }
  
  beta_hat <- sol$x
  
  SA_pi_hat <- plogis(as.numeric(SA_z %*% beta_hat))
  SB_pi_hat <- plogis(as.numeric(SB_z %*% beta_hat))
  
  N_hat_A <- sum(1 / SA_pi_hat)
  N_hat_B <- sum(SB_d)
  
  A_hat <- t(SA_z) %*% ((1 - SA_pi_hat) * SA_z) / N_hat_A
  
  A_mu_hat <- colSums(
    SA_z * as.numeric((1 - SA_pi_hat) * SA_mu)
  ) / N_hat_A
  
  A_inv_A_mu <- tryCatch(
    qr.solve(A_hat, A_mu_hat),
    error = function(e) rep(NA_real_, length(A_mu_hat))
  )
  
  if (any(is.na(A_inv_A_mu))) {
    next
  }
  
  part_A <- as.numeric(crossprod(A_mu_hat, A_inv_A_mu))
  
  term_beta <- as.numeric(SB_z %*% A_inv_A_mu)
  
  # ============================================================
  # Construction (a): S1
  # ============================================================
  
  A2_hat_s1 <- sum((1 - SA_pi_hat)^2 * (SA_mu^2 + SA_sigma2)) / N_hat_A
  
  SA_xx <- SA_x[, 2:3, drop = FALSE]
  
  a_hat <- t(as.numeric(1 / SA_sigma2) * SA_xx) %*% SA_xx / N_hat_A
  
  b_hat <- matrix(0.5 * nrow(SA_xx) / N_hat_A, nrow = 1, ncol = 1)
  
  B_hat <- as.matrix(Matrix::bdiag(a_hat, b_hat))
  
  B1_hat <- colSums(
    cbind(SA_xx, rep(0, length(SA_y))) * as.numeric(1 - SA_pi_hat)
  ) / N_hat_A
  
  part_B <- tryCatch(
    as.numeric(crossprod(B1_hat, qr.solve(B_hat, B1_hat))),
    error = function(e) NA_real_
  )
  
  B2_hat <- sum(SA_pi_hat * (1 - SA_pi_hat) * SA_mu^2) / N_hat_A
  
  sigma1_s1 <- A2_hat_s1 + B2_hat - part_A - part_B
  
  phi_hat_s1 <- SB_pi_hat * (SB_mu - term_beta)
  
  sigma_HT_s1 <- sum(SB_d * (SB_d - 1) * phi_hat_s1^2) / N_hat_B
  
  sigma_hat_s1 <- sigma1_s1 + sigma_HT_s1
  
  if (is.na(sigma1_s1) || sigma1_s1 <= 0 || 
      is.na(sigma_hat_s1) || sigma_hat_s1 <= 0) {
    
    count_neg_s1 <- count_neg_s1 + 1
    
  } else {
    
    S1_stat <- sum(SA_y * (1 - SA_pi_hat) + SA_pi_hat * SA_mu) -
      sum(SB_d * SB_pi_hat * SB_mu)
    
    Z1_stat <- abs(S1_stat) / sqrt(N * sigma_hat_s1)
    
    p_values_s1[b] <- 2 * (1 - pnorm(Z1_stat))
  }
  
  # ============================================================
  # Construction (b): S2
  # ============================================================
  
  A4_hat_s2 <- sum((1 - SA_pi_hat) * SA_mu^2) / N_hat_A
  
  sigma3_s2 <- A4_hat_s2 - part_A
  
  phi_hat_s2 <- SB_pi_hat * (SB_mu - term_beta)
  
  sigma_HT_s2 <- sum(SB_d * (SB_d - 1) * phi_hat_s2^2) / N_hat_B
  
  sigma_hat_s2 <- sigma3_s2 + sigma_HT_s2
  
  if (is.na(sigma3_s2) || sigma3_s2 <= 0 || 
      is.na(sigma_hat_s2) || sigma_hat_s2 <= 0) {
    
    count_neg_s2 <- count_neg_s2 + 1
    
  } else {
    
    S2_stat <- sum(SA_mu) - sum(SB_d * SB_pi_hat * SB_mu)
    
    Z2_stat <- abs(S2_stat) / sqrt(N * sigma_hat_s2)
    
    p_values_s2[b] <- 2 * (1 - pnorm(Z2_stat))
  }
}

cat("\nSimulation finished.\n")
cat(sprintf("beta fail count: %d\n", count_fail_beta))
cat(sprintf("S1 sigma <= 0 count: %d\n", count_neg_s1))
cat(sprintf("S2 sigma <= 0 count: %d\n", count_neg_s2))
cat(sprintf("valid S1 p-values: %d\n", sum(!is.na(p_values_s1))))
cat(sprintf("valid S2 p-values: %d\n", sum(!is.na(p_values_s2))))

# ============================================================
# Size calibration over alpha = 0.01,...,0.30
# ============================================================

alpha_grid <- seq(0.01, 0.30, by = 0.01)

size_s1 <- sapply(alpha_grid, function(a) {
  mean(p_values_s1 <= a, na.rm = TRUE)
})

size_s2 <- sapply(alpha_grid, function(a) {
  mean(p_values_s2 <= a, na.rm = TRUE)
})

size_results <- data.frame(
  alpha = alpha_grid,
  S1_Size = size_s1,
  S2_Size = size_s2
)

print(size_results)

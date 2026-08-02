rm(list = ls())
library(nleqslv)
library(Matrix)
library(ggplot2)

set.seed(2026)

N <- 100000
nA <- 2000
nB <- 2000
sim_times <- 2000

# ============================================================
# Case A: Y | X ~ N(-X + X^2, exp(0.5))
# ============================================================

x <- rnorm(N, mean = 1, sd = 1)

X_pop <- cbind(1, x, x^2)   # outcome model: (1, x, x^2)
Z_pop <- cbind(1, x)        # inclusion model: (1, x)

lin_pred <- -1 * x + 1 * x^2
lambda_true <- exp(0.5)

Y <- lin_pred + sqrt(lambda_true) * rnorm(N, 0, 1)

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
pi_B_pop <- pmin(pi_B_pop, 1)

# ============================================================
# Outcome model MLE
# Y | X ~ N(xi1*x + xi2*x^2, exp(xi3))
# ============================================================

xiloglik <- function(xi, x, y, d) {
  sum(d * (xi[3] + (y - xi[1] * x - xi[2] * x^2)^2 * exp(-xi[3])))
}

max.xiloglik <- function(xi, x, y, d) {
  nlminb(xi, xiloglik, x = x, y = y, d = d)$par
}

# ============================================================
# H0: gamma = 0
# ============================================================

gamma_val <- 0

find_theta0 <- function(be0) {
  sum(plogis(be0 + 0.5 * x + gamma_val * Y)) - nA
}

beta0_opt <- uniroot(find_theta0, c(-20, 20))$root

pi_A_pop <- plogis(beta0_opt + 0.5 * x + gamma_val * Y)

cat(sprintf("H0: gamma = %.2f\n", gamma_val))
cat(sprintf("beta0 = %.4f, sum(pi_A) = %.2f\n", beta0_opt, sum(pi_A_pop)))

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
  # Estimate outcome model
  # ------------------------------------------------------------
  
  xi_init <- c(-1, 1, 0.5)
  
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
    sum_SB <- colSums(SB_d * pi_B_est * SB_z)
    
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
  
  # Common A matrix
  A_hat <- t(SA_z) %*% ((1 - SA_pi_hat) * SA_z) / N_hat_A
  
  # ============================================================
  # Construction (a): S1
  # ============================================================
  
  A1_hat_s1 <- t(SA_z) %*% ((1 - SA_pi_hat) * SA_mu) / N_hat_A
  A2_hat_s1 <- sum((1 - SA_pi_hat)^2 * (SA_mu^2 + SA_sigma2)) / N_hat_A
  
  SA_xx <- SA_x[, 2:3, drop = FALSE]
  
  a_hat <- t(as.numeric(1 / SA_sigma2) * SA_xx) %*% SA_xx / N_hat_A
  
  b_hat <- matrix(0.5 * nrow(SA_xx) / N_hat_A, nrow = 1, ncol = 1)
  
  B_hat <- as.matrix(Matrix::bdiag(a_hat, b_hat))
  
  B1_hat_s1 <- t(as.numeric(1 - SA_pi_hat)) %*%
    cbind(SA_xx, rep(0, length(SA_y))) / N_hat_A
  
  B2_hat_s1 <- sum(SA_pi_hat * (1 - SA_pi_hat) * SA_mu^2) / N_hat_A
  
  part_A_s1 <- tryCatch(
    as.numeric(t(A1_hat_s1) %*% solve(A_hat) %*% A1_hat_s1),
    error = function(e) NA_real_
  )
  
  part_B_s1 <- tryCatch(
    as.numeric(B1_hat_s1 %*% solve(B_hat) %*% t(B1_hat_s1)),
    error = function(e) NA_real_
  )
  
  sigma1_s1 <- A2_hat_s1 + B2_hat_s1 - part_A_s1 - part_B_s1
  
  if (is.na(sigma1_s1) || sigma1_s1 <= 0) {
    count_neg_s1 <- count_neg_s1 + 1
  } else {
    
    term2_s1 <- as.numeric(SB_z %*% solve(A_hat) %*% A1_hat_s1)
    
    phi_hat_s1 <- SB_pi_hat * SB_mu - SB_pi_hat * term2_s1
    
    sigma_HT_s1 <- sum(SB_d * (SB_d - 1) * phi_hat_s1^2) / N_hat_B
    
    sigma_hat_s1 <- sigma1_s1 + sigma_HT_s1
    
    if (is.na(sigma_hat_s1) || sigma_hat_s1 <= 0) {
      count_neg_s1 <- count_neg_s1 + 1
    } else {
      
      S1_stat <- sum(SA_y * (1 - SA_pi_hat) + SA_pi_hat * SA_mu) -
        sum(SB_d * SB_pi_hat * SB_mu)
      
      Z1_stat <- abs(S1_stat) / sqrt(N * sigma_hat_s1)
      
      p_values_s1[b] <- 2 * (1 - pnorm(Z1_stat))
    }
  }
  
  # ============================================================
  # Construction (b): S2
  # ============================================================
  
  A3_hat_s2 <- t(SA_z) %*% ((1 - SA_pi_hat) * SA_mu) / N_hat_A
  A4_hat_s2 <- sum((1 - SA_pi_hat) * SA_mu^2) / N_hat_A
  
  part_A_s2 <- tryCatch(
    as.numeric(t(A3_hat_s2) %*% solve(A_hat) %*% A3_hat_s2),
    error = function(e) NA_real_
  )
  
  sigma3_s2 <- A4_hat_s2 - part_A_s2
  
  if (is.na(sigma3_s2) || sigma3_s2 <= 0) {
    count_neg_s2 <- count_neg_s2 + 1
  } else {
    
    term2_s2 <- as.numeric(SB_z %*% solve(A_hat) %*% A3_hat_s2)
    
    phi_hat_s2 <- SB_pi_hat * (term2_s2 - SB_mu)
    
    sigma_HT_s2 <- sum(SB_d * (SB_d - 1) * phi_hat_s2^2) / N_hat_B
    
    sigma_hat_s2 <- sigma3_s2 + sigma_HT_s2
    
    if (is.na(sigma_hat_s2) || sigma_hat_s2 <= 0) {
      count_neg_s2 <- count_neg_s2 + 1
    } else {
      
      S2_stat <- sum(SA_mu) - sum(SB_d * SB_pi_hat * SB_mu)
      
      Z2_stat <- abs(S2_stat) / sqrt(N * sigma_hat_s2)
      
      p_values_s2[b] <- 2 * (1 - pnorm(Z2_stat))
    }
  }
}

cat("\nSimulation finished.\n")
cat(sprintf("beta fail count: %d\n", count_fail_beta))
cat(sprintf("S1 sigma <= 0 count: %d\n", count_neg_s1))
cat(sprintf("S2 sigma <= 0 count: %d\n", count_neg_s2))

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

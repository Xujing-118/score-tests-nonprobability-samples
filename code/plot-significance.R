rm(list = ls())

library(nleqslv)
library(Matrix)
library(ggplot2)
setwd("D:/paper/MAR/simulation-mis")

lambda <- 0
#lambda <- 0.1

N <- 100000
n_grid <- c(500, 1000, 2000)
sim_times <- 2000
gamma_val <- 0
alpha_grid <- seq(0.01, 0.30, by = 0.01)

stopifnot(is.numeric(lambda), length(lambda) == 1, lambda >= 0)

lambda_tag <- if (lambda == 0) "0" else "0.1"
output_dir <- file.path(getwd(), paste0("figuresin_lambda_", lambda_tag))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Y | X ~ N(xi1*x + xi2*x^2, exp(xi3))
xiloglik <- function(xi, x, y, d) {
  sum(d * (xi[3] + (y - xi[1] * x - xi[2] * x^2)^2 * exp(-xi[3])))
}

max.xiloglik <- function(xi, x, y, d) {
  nlminb(xi, xiloglik, x = x, y = y, d = d)$par
}

# Case A : Y | X ~ N(-X + X^2, exp(0.5))

run_case_A <- function(nA, nB) {
  set.seed(2026)

  x <- rnorm(N, mean = 0, sd = 1)
  X_pop <- cbind(1, x, x^2)

  lin_pred <- -1 * x + 1 * x^2
  lambda_true <- exp(0.5)

  Y <- lin_pred + sqrt(lambda_true) * rnorm(N, 0, 1)

  temp_val <- x + 0.03 * Y
  c_val <- (max(temp_val) - 50 * min(temp_val)) / 49

  if (c_val + min(temp_val) <= 0) {
    c_val <- abs(min(temp_val)) + 0.1
  }

  z_size <- c_val + temp_val
  pi_B_pop <- nB * z_size / sum(z_size)
  pi_B_pop <- pmin(pi_B_pop, 1)

  # H0: gamma = 0

  find_theta0 <- function(be0) {
    sum(plogis(be0 + 0.5 * x + lambda * sin(x) + gamma_val * Y)) - nA
  }

  beta0_opt <- uniroot(find_theta0, c(-20, 20))$root

  pi_A_pop <- plogis(beta0_opt + 0.5 * x + lambda * sin(x) + gamma_val * Y)

  cat(sprintf(
    "\nCase A: nA=nB=%d, lambda=%.1f, beta0=%.4f, sum(pi_A)=%.2f\n",
    nA, lambda, beta0_opt, sum(pi_A_pop)
  ))

  set.seed(2026)
  
  p_values_s1 <- rep(NA_real_, sim_times)
  p_values_s2 <- rep(NA_real_, sim_times)

  count_neg_s1 <- 0
  count_neg_s2 <- 0
  count_fail_beta <- 0

  for (b in 1:sim_times) {

    if (b %% 500 == 0) {
      cat(sprintf("  %d / %d ...\n", b, sim_times))
    }

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

    # Outcome model
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

    A_hat <- t(SA_z) %*% ((1 - SA_pi_hat) * SA_z) / N_hat_A
    A1_hat_s1 <- t(SA_z) %*% ((1 - SA_pi_hat) * SA_mu) / N_hat_A
    A2_hat_s1 <- sum( (1 - SA_pi_hat)^2 * (SA_mu^2 + SA_sigma2) ) / N_hat_A

    SA_xx <- SA_x[, 2:3, drop = FALSE]
    a_hat <- t(as.numeric(1 / SA_sigma2) * SA_xx) %*% SA_xx / N_hat_A
    b_hat <- matrix(
      0.5 * nrow(SA_xx) / N_hat_A,
      nrow = 1,
      ncol = 1
    )

    B_hat <- as.matrix(Matrix::bdiag(a_hat, b_hat))

    B1_hat_s1 <- t(as.numeric(1 - SA_pi_hat)) %*%cbind(SA_xx, rep(0, length(SA_y))) / N_hat_A
    B2_hat_s1 <- sum( SA_pi_hat * (1 - SA_pi_hat) * SA_mu^2 ) / N_hat_A

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
      term2_s1 <- as.numeric( SB_z %*% solve(A_hat) %*% A1_hat_s1 )

      phi_hat_s1 <- SB_pi_hat * SB_mu - SB_pi_hat * term2_s1

      sigma_HT_s1 <- sum(SB_d * (SB_d - 1) * phi_hat_s1^2) / N_hat_B
      sigma_hat_s1 <- sigma1_s1 + sigma_HT_s1

      if (is.na(sigma_hat_s1) || sigma_hat_s1 <= 0) {
        count_neg_s1 <- count_neg_s1 + 1
      } else {
        S1_stat <- sum( SA_y * (1 - SA_pi_hat) + SA_pi_hat * SA_mu ) - sum(SB_d * SB_pi_hat * SB_mu)
        Z1_stat <- abs(S1_stat) / sqrt(N * sigma_hat_s1)
        p_values_s1[b] <- 2 * (1 - pnorm(Z1_stat))
      }
    }

    A3_hat_s2 <- t(SA_z) %*%((1 - SA_pi_hat) * SA_mu) / N_hat_A
    A4_hat_s2 <- sum( (1 - SA_pi_hat) * SA_mu^2 ) / N_hat_A

    part_A_s2 <- tryCatch(
      as.numeric(t(A3_hat_s2) %*% solve(A_hat) %*% A3_hat_s2),
      error = function(e) NA_real_
    )

    sigma3_s2 <- A4_hat_s2 - part_A_s2

    if (is.na(sigma3_s2) || sigma3_s2 <= 0) {
      count_neg_s2 <- count_neg_s2 + 1
    } else {

      term2_s2 <- as.numeric( SB_z %*% solve(A_hat) %*% A3_hat_s2)
      phi_hat_s2 <- SB_pi_hat * (term2_s2 - SB_mu)

      sigma_HT_s2 <- sum(SB_d * (SB_d - 1) * phi_hat_s2^2 ) / N_hat_B
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

  cat(sprintf(
    "Case A n=%d finished: beta fail=%d, S1 invalid=%d, S2 invalid=%d\n",
    nA, count_fail_beta, count_neg_s1, count_neg_s2
  ))

  size_s1 <- sapply(alpha_grid, function(a) {
    mean(p_values_s1 <= a, na.rm = TRUE)
  })

  size_s2 <- sapply(alpha_grid, function(a) {
    mean(p_values_s2 <= a, na.rm = TRUE)
  })

  list(
    size = data.frame(
      Case = "A",
      n = nA,
      alpha = alpha_grid,
      S1_Size = size_s1,
      S2_Size = size_s2
    ),
    p = data.frame(
      Case = "A",
      n = nA,
      replicate = 1:sim_times,
      p_value_s1 = p_values_s1,
      p_value_s2 = p_values_s2
    )
  )
}


# Case B
# Y | X, R=1 ~ N(-X + X^2, exp(2))

run_case_B <- function(nA, nB) {

  set.seed(2026)

  x <- rnorm(N, mean = 0, sd = 1)
  X_pop <- cbind(1, x, x^2)

  mu_A <- -1 * x + 1 * x^2
  sigma2_true <- exp(2)
  sigma_true <- sqrt(sigma2_true)

  # ----------------------------------------------------------
  # H0: gamma = 0
  # 保持原 Case B 构造，只在 participation linear predictor
  # 中增加 lambda*x^2。
  # ----------------------------------------------------------

  c_gamma <- -gamma_val * mu_A + 0.5 * gamma_val^2 * sigma2_true

  beta0_opt <- uniroot(function(be0) {
    sum(plogis( be0 + 0.5 * x + lambda * sin(x) - c_gamma ) ) - nA
  }, c(-20, 20))$root

  pi_x <- plogis( beta0_opt + 0.5 * x + lambda * sin(x) - c_gamma)

  set.seed(2026)

  G <- rbinom(N, 1, pi_x)
  Y <- rnorm(
    N,
    mean = ifelse(
      G == 1,
      mu_A,
      mu_A - gamma_val * sigma2_true
    ),
    sd = sigma_true
  )

  pi_A_pop <- plogis( beta0_opt + 0.5 * x + lambda * sin(x) + gamma_val * Y )

  cat(sprintf(
    "\nCase B: nA=nB=%d, lambda=%.1f, beta0=%.4f, sum(pi_A)=%.2f\n",
    nA, lambda, beta0_opt, sum(pi_A_pop)
  ))

  temp_val <- x + 0.03 * Y
  c_val <- (max(temp_val) - 50 * min(temp_val)) / 49

  if (c_val + min(temp_val) <= 0) {
    c_val <- abs(min(temp_val)) + 0.1
  }

  z_size <- c_val + temp_val
  pi_B_pop <- nB * z_size / sum(z_size)
  pi_B_pop <- pmin(pi_B_pop, 0.99)

  p_values_s1 <- rep(NA_real_, sim_times)
  p_values_s2 <- rep(NA_real_, sim_times)

  count_neg_s1 <- 0
  count_neg_s2 <- 0
  count_fail_beta <- 0

  for (b in 1:sim_times) {

    if (b %% 500 == 0) {
      cat(sprintf("  %d / %d ...\n", b, sim_times))
    }

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

    xi_init <- c(-1, 1, 2)
    xi_hat <- max.xiloglik(
      xi_init,
      SA_x[, 2],
      SA_y,
      rep(1, length(SA_y))
    )

    SA_mu <- as.numeric( xi_hat[1] * SA_x[, 2] +  xi_hat[2] * SA_x[, 3] )
    SB_mu <- as.numeric( xi_hat[1] * SB_x[, 2] + xi_hat[2] * SB_x[, 3] )

    SA_sigma2 <- exp(xi_hat[3])

    SA_z <- SA_x[, 1:2, drop = FALSE]
    SB_z <- SB_x[, 1:2, drop = FALSE]

    score_eqn <- function(beta) {
      pi_B_est <- plogis(as.numeric(SB_z %*% beta))

      sum_SA <- colSums(SA_z)
      sum_SB <- colSums(SB_z * as.numeric(SB_d * pi_B_est) )
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

    A_hat <- t(SA_z) %*%((1 - SA_pi_hat) * SA_z) / N_hat_A
    A_mu_hat <- colSums( SA_z * as.numeric((1 - SA_pi_hat) * SA_mu)) / N_hat_A

    A_inv_A_mu <- tryCatch(
      qr.solve(A_hat, A_mu_hat),
      error = function(e) rep(NA_real_, length(A_mu_hat))
    )

    if (any(is.na(A_inv_A_mu))) {
      next
    }

    part_A <- as.numeric(crossprod(A_mu_hat, A_inv_A_mu))
    term_beta <- as.numeric(SB_z %*% A_inv_A_mu)

    A2_hat_s1 <- sum( (1 - SA_pi_hat)^2 *  (SA_mu^2 + SA_sigma2) ) / N_hat_A
    SA_xx <- SA_x[, 2:3, drop = FALSE]

    a_hat <- t(as.numeric(1 / SA_sigma2) * SA_xx) %*% SA_xx / N_hat_A

    b_hat <- matrix(
      0.5 * nrow(SA_xx) / N_hat_A,
      nrow = 1,
      ncol = 1
    )

    B_hat <- as.matrix(Matrix::bdiag(a_hat, b_hat))

    B1_hat <- colSums(cbind(SA_xx, rep(0, length(SA_y))) * as.numeric(1 - SA_pi_hat)) / N_hat_A

    part_B <- tryCatch(
      as.numeric(
        crossprod(
          B1_hat,
          qr.solve(B_hat, B1_hat)
        )
      ),
      error = function(e) NA_real_
    )

    B2_hat <- sum(SA_pi_hat * (1 - SA_pi_hat) * SA_mu^2) / N_hat_A

    sigma1_s1 <- A2_hat_s1 + B2_hat - part_A - part_B
    phi_hat_s1 <- SB_pi_hat * (SB_mu - term_beta)

    sigma_HT_s1 <- sum( SB_d * (SB_d - 1) * phi_hat_s1^2) / N_hat_B
    sigma_hat_s1 <- sigma1_s1 + sigma_HT_s1

    if (
      is.na(sigma1_s1) || sigma1_s1 <= 0 ||
      is.na(sigma_hat_s1) || sigma_hat_s1 <= 0
    ) {

      count_neg_s1 <- count_neg_s1 + 1

    } else {

      S1_stat <- sum(  SA_y * (1 - SA_pi_hat) + SA_pi_hat * SA_mu ) - sum(SB_d * SB_pi_hat * SB_mu)
      Z1_stat <- abs(S1_stat) / sqrt(N * sigma_hat_s1)
      p_values_s1[b] <- 2 * (1 - pnorm(Z1_stat))
    }

    A4_hat_s2 <- sum( (1 - SA_pi_hat) * SA_mu^2 ) / N_hat_A
    sigma3_s2 <- A4_hat_s2 - part_A

    phi_hat_s2 <- SB_pi_hat * (SB_mu - term_beta)

    sigma_HT_s2 <- sum( SB_d * (SB_d - 1) * phi_hat_s2^2 ) / N_hat_B
    sigma_hat_s2 <- sigma3_s2 + sigma_HT_s2

    if (
      is.na(sigma3_s2) || sigma3_s2 <= 0 ||
      is.na(sigma_hat_s2) || sigma_hat_s2 <= 0
    ) {
      count_neg_s2 <- count_neg_s2 + 1
    } else {

      S2_stat <- sum(SA_mu) - sum(SB_d * SB_pi_hat * SB_mu)
      Z2_stat <- abs(S2_stat) / sqrt(N * sigma_hat_s2)

      p_values_s2[b] <- 2 * (1 - pnorm(Z2_stat))
    }
  }

  cat(sprintf(
    "Case B n=%d finished: beta fail=%d, S1 invalid=%d, S2 invalid=%d\n",
    nA, count_fail_beta, count_neg_s1, count_neg_s2
  ))

  size_s1 <- sapply(alpha_grid, function(a) {
    mean(p_values_s1 <= a, na.rm = TRUE)
  })

  size_s2 <- sapply(alpha_grid, function(a) {
    mean(p_values_s2 <= a, na.rm = TRUE)
  })

  list(
    size = data.frame(
      Case = "B",
      n = nA,
      alpha = alpha_grid,
      S1_Size = size_s1,
      S2_Size = size_s2
    ),
    p = data.frame(
      Case = "B",
      n = nA,
      replicate = 1:sim_times,
      p_value_s1 = p_values_s1,
      p_value_s2 = p_values_s2
    )
  )
}

results <- list()
kk <- 0

for (n in n_grid) {
  kk <- kk + 1
  results[[kk]] <- run_case_A(nA = n, nB = n)
}

for (n in n_grid) {
  kk <- kk + 1
  results[[kk]] <- run_case_B(nA = n, nB = n)
}

size_all <- do.call(
  rbind,
  lapply(results, function(z) z$size)
)

p_all <- do.call(
  rbind,
  lapply(results, function(z) z$p)
)

write.csv(
  size_all,
  file.path(output_dir, "size_all.csv"),
  row.names = FALSE
)

write.csv(
  p_all,
  file.path(output_dir, "p_all.csv"),
  row.names = FALSE
)


# Plot

plot_df <- rbind(
  data.frame(
    alpha = size_all$alpha,
    Size = size_all$S1_Size,
    Method = "Z1",
    Panel = paste0(size_all$Case, " with n=", size_all$n)
  ),
  data.frame(
    alpha = size_all$alpha,
    Size = size_all$S2_Size,
    Method = "Z2",
    Panel = paste0(size_all$Case, " with n=", size_all$n)
  )
)

plot_df$Method <- factor(
  plot_df$Method,
  levels = c("Z1", "Z2")
)

plot_df$Panel <- factor(
  plot_df$Panel,
  levels = c(
    "A with n=500",
    "A with n=1000",
    "A with n=2000",
    "B with n=500",
    "B with n=1000",
    "B with n=2000"
  )
)

p <- ggplot(
  plot_df,
  aes(
    x = alpha,
    y = Size,
    color = Method,
    shape = Method,
    group = Method
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = "gray60",
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  # geom_line(linewidth = 0.8, linetype = "solid") +
  geom_point(size = 2, stroke = 0.9) +
  facet_wrap(~ Panel, ncol = 3) +
  scale_color_manual(
    values = c(
      "Z1" = "#B22222",
      "Z2" = "#1F4E79"
    ),
    labels = c(
      "Z1" = expression(Z[1]),
      "Z2" = expression(Z[2])
    )
  ) +
  scale_shape_manual(
    values = c(
      "Z1" = 16,
      "Z2" = 3
    ),
    labels = c(
      "Z1" = expression(Z[1]),
      "Z2" = expression(Z[2])
    )
  ) +
  scale_x_continuous(
    breaks = seq(0, 0.30, by = 0.05),
    limits = c(0, 0.33),
    labels = function(x) sprintf("%.2f", x)
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.30, by = 0.05),
    limits = c(0, 0.33),
    labels = function(x) sprintf("%.2f", x)
  ) +
  coord_equal() +
  labs(
    x = "Significance level",
    y = "Empirical rejection rate"
  ) +
  theme_bw(base_size = 15) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 16,
      color = "black"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "black",
      linewidth = 0.65
    ),
    axis.text = element_text(
      color = "black",
      size = 11
    ),
    axis.title = element_text(
      color = "black",
      size = 17
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.6
    ),
    axis.ticks.length = grid::unit(0.18, "cm"),

    legend.position = "inside",
    legend.position.inside = c(0.005, 0.98),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(
      size = 13,
      color = "black"
    ),
    legend.key.width = grid::unit(0.7, "cm"),
    legend.key.height = grid::unit(0.45, "cm"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),

    plot.margin = margin(8, 8, 8, 8),
    panel.spacing = grid::unit(1.1, "lines")
  )

print(p)

figure_name <- if (lambda == 0) {
  "Figuresin1_lambda_0.png"
} else {
  "Figuresin2_lambda_0.1.png"
}

ggsave(
  filename = file.path(output_dir, figure_name),
  plot = p,
  width = 10.5,
  height = 7.2,
  dpi = 600
)

cat(sprintf(
  "\nFinished. lambda=%.1f\nResults saved in: %s\n",
  lambda,
  normalizePath(output_dir, winslash = "/")
))

########################
##########################
############################
###########################

suppressWarnings({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
})
library(ggplot2)

# ============================================================
# One-pass averaged SGD + perturbed SGD for quantile regression
# ============================================================
run_qr_sgd_once <- function(
    X,
    Y,
    B,
    tau = 0.5,
    step_c = 1,
    step_alpha = 2 / 3,
    burn_in = 0L
) {
  N <- nrow(X)
  p <- ncol(X)

  theta_hat <- numeric(p)
  theta_bar <- numeric(p)

  theta_star_hat <- matrix(0, nrow = B, ncol = p)
  theta_star_bar <- matrix(0, nrow = B, ncol = p)

  n_keep <- 0L

  for (n in seq_len(N)) {
    x_n <- X[n, ]
    y_n <- Y[n]
    gamma_n <- step_c * n^(-step_alpha)

    # Main SGD update
    resid_main <- y_n - sum(x_n * theta_hat)
    score_main <- tau - as.numeric(resid_main < 0)
    theta_hat <- theta_hat + gamma_n * score_main * x_n

    # Perturbed SGD updates
    resid_star <- y_n - as.vector(theta_star_hat %*% x_n)
    score_star <- tau - as.numeric(resid_star < 0)
    w_n <- rexp(B, rate = 1)

    theta_star_hat <- theta_star_hat +
      gamma_n * ((w_n * score_star) %o% x_n)

    # Tail averaging after burn-in
    if (n > burn_in) {
      n_keep <- n_keep + 1L
      theta_bar <- theta_bar + (theta_hat - theta_bar) / n_keep
      theta_star_bar <- theta_star_bar +
        (theta_star_hat - theta_star_bar) / n_keep
    }
  }

  list(
    theta_bar = theta_bar,
    theta_star_bar = theta_star_bar
  )
}

# ============================================================
# CI constructors for the first coordinate
# ============================================================
usual_sgd_ci <- function(theta_bar_1, theta_star_bar_1, alpha = 0.10) {
  B <- length(theta_star_bar_1)
  z <- sort(theta_star_bar_1)

  left_idx <- ceiling(B * alpha / 2)
  right_idx <- ceiling(B * (1 - alpha / 2))

  left_idx <- max(1L, min(B, left_idx))
  right_idx <- max(1L, min(B, right_idx))

  L <- 2 * theta_bar_1 - z[right_idx]
  U <- 2 * theta_bar_1 - z[left_idx]

  c(L = L, U = U)
}

modified_sgd_ci <- function(theta_bar_1, theta_star_bar_1, alpha = 0.10) {
  B <- length(theta_star_bar_1)
  z <- sort(theta_star_bar_1)

  left_idx <- floor((B + 1) * alpha / 2)
  right_idx <- ceiling((1 - alpha) * (B + 1)) + floor((alpha / 2) * (B + 1))

  left_idx <- max(1L, min(B, left_idx))
  right_idx <- max(1L, min(B, right_idx))

  L <- 2 * theta_bar_1 - z[right_idx]
  U <- 2 * theta_bar_1 - z[left_idx]

  c(L = L, U = U)
}

# ============================================================
# Main Monte Carlo comparison (parallelized over n_iter)
# ============================================================
simulate_qr_sgd_ci <- function(
    B_grid = seq(1, 200),
    n_iter = 1000,
    N = 10000,
    alpha = 0.10,
    tau = 0.50,
    theta0 = c(0.2, -0.2, 0),
    step_c = 1,
    step_alpha = 2 / 3,
    burn_in = NULL,
    seed = 123,
    n_cores = parallel::detectCores(logical = TRUE)
) {
  set.seed(seed)

  stopifnot(length(theta0) == 3L)

  if (is.null(burn_in)) {
    if (N == 10000) {
      burn_in <- 2000L
    } else if (N == 20000) {
      burn_in <- 4000L
    } else {
      burn_in <- floor(0.2 * N)
    }
  }

  n_cores <- as.integer(n_cores)
  if (is.na(n_cores) || n_cores < 1L) {
    n_cores <- 1L
  }

  B_max <- max(B_grid)
  theta01_true <- theta0[1]
  n_B <- length(B_grid)

  iter_results <- parallel::mclapply(
    X = seq_len(n_iter),
    FUN = function(iter) {
      X <- matrix(rnorm(N * 3), nrow = N, ncol = 3)
      eps <- rexp(N, rate = 1) - rexp(N, rate = 1)
      Y <- as.vector(X %*% theta0 + eps)

      fit <- run_qr_sgd_once(
        X = X,
        Y = Y,
        B = B_max,
        tau = tau,
        step_c = step_c,
        step_alpha = step_alpha,
        burn_in = burn_in
      )

      theta_bar_1 <- fit$theta_bar[1]
      theta_star_all_1 <- fit$theta_star_bar[, 1]

      covered_usual <- integer(n_B)
      covered_modified <- integer(n_B)
      widths_usual <- numeric(n_B)
      widths_modified <- numeric(n_B)

      for (j in seq_along(B_grid)) {
        B <- B_grid[j]
        theta_star_1 <- theta_star_all_1[seq_len(B)]

        ci_usual <- usual_sgd_ci(
          theta_bar_1 = theta_bar_1,
          theta_star_bar_1 = theta_star_1,
          alpha = alpha
        )

        ci_modified <- modified_sgd_ci(
          theta_bar_1 = theta_bar_1,
          theta_star_bar_1 = theta_star_1,
          alpha = alpha
        )

        widths_usual[j] <- ci_usual["U"] - ci_usual["L"]
        widths_modified[j] <- ci_modified["U"] - ci_modified["L"]

        covered_usual[j] <- as.integer(theta01_true >= ci_usual["L"] &&
                                         theta01_true < ci_usual["U"])
        covered_modified[j] <- as.integer(theta01_true >= ci_modified["L"] &&
                                            theta01_true < ci_modified["U"])
      }

      list(
        covered_usual = covered_usual,
        covered_modified = covered_modified,
        widths_usual = widths_usual,
        widths_modified = widths_modified
      )
    },
    mc.cores = n_cores
  )

  covered_usual <- Reduce(
    f = `+`,
    x = lapply(iter_results, function(x) x$covered_usual)
  )
  covered_modified <- Reduce(
    f = `+`,
    x = lapply(iter_results, function(x) x$covered_modified)
  )

  widths_usual <- do.call(
    rbind,
    lapply(iter_results, function(x) x$widths_usual)
  )
  widths_modified <- do.call(
    rbind,
    lapply(iter_results, function(x) x$widths_modified)
  )

  out <- data.frame(
    B = integer(),
    method = character(),
    coverage = numeric(),
    mc_se = numeric(),
    mean_width = numeric(),
    width_se = numeric(),
    stringsAsFactors = FALSE
  )

  for (j in seq_along(B_grid)) {
    B <- B_grid[j]

    cov_u <- covered_usual[j] / n_iter
    cov_m <- covered_modified[j] / n_iter

    se_u <- sqrt(cov_u * (1 - cov_u) / n_iter)
    se_m <- sqrt(cov_m * (1 - cov_m) / n_iter)

    mw_u <- mean(widths_usual[, j])
    mw_m <- mean(widths_modified[, j])

    wse_u <- sd(widths_usual[, j]) / sqrt(n_iter)
    wse_m <- sd(widths_modified[, j]) / sqrt(n_iter)

    out <- rbind(
      out,
      data.frame(
        B = B,
        method = "Usual SGD CI",
        coverage = cov_u,
        mc_se = se_u,
        mean_width = mw_u,
        width_se = wse_u
      ),
      data.frame(
        B = B,
        method = "Modified SGD CI",
        coverage = cov_m,
        mc_se = se_m,
        mean_width = mw_m,
        width_se = wse_m
      )
    )
  }

  out
}

# ============================================================
# Plot helpers
# ============================================================
.make_ticks <- function(B_vals, n_ticks = 10) {
  unique(as.integer(seq(min(B_vals), max(B_vals), length.out = n_ticks)))
}

plot_coverage <- function(df, alpha, n_ticks = 10) {
  tick_breaks <- .make_ticks(df$B, n_ticks)

  ggplot(df, aes(x = B, y = coverage, color = method, group = method)) +
    geom_hline(
      yintercept = 1 - alpha,
      linetype = "dotted",
      color = "black",
      linewidth = 1.2
    ) +
    geom_line() +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = tick_breaks) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(
      title = expression(paste("Coverage for CI of ", theta[0, 1])),
      y = "Estimated Coverage",
      x = "Number of perturbation replicates (B)",
      color = "Method"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 13),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5)
    )
}

plot_width <- function(df, n_ticks = 10) {
  tick_breaks <- .make_ticks(df$B, n_ticks)

  ggplot(df, aes(x = B, y = mean_width, color = method, group = method)) +
    geom_line() +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = tick_breaks) +
    labs(
      title = expression(paste("Mean CI width for ", theta[0, 1])),
      y = "Mean Width",
      x = "Number of perturbation replicates (B)",
      color = "Method"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 13),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5)
    )
}

# ============================================================
# Example usage
# ============================================================
res_df <- simulate_qr_sgd_ci(
  B_grid = seq(1, 200),
  n_iter = 1000,
  N = 10000,
  alpha = 0.10,
  tau = 0.50,
  theta0 = c(0.2, -0.2, 0),
  step_c = 1,
  step_alpha = 2 / 3,
  seed = 123
)

save(res_df, file = "qr_sgd_ci_compare.RData")

p_cov <- plot_coverage(res_df, alpha = 0.10, n_ticks = 10)
p_wid <- plot_width(res_df, n_ticks = 10)

dir.create("PDF_images", recursive = TRUE, showWarnings = FALSE)
ggsave(filename = file.path("PDF_images", "sgd_cov.pdf"), plot = p_cov, device = "pdf", width = 6, height = 4)
ggsave(filename = file.path("PDF_images", "sgd_width.pdf"), plot = p_wid, device = "pdf", width = 6, height = 4)

print(p_cov)
print(p_wid)

head(res_df)

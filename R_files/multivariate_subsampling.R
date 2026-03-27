suppressWarnings({
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
})
library(ggplot2)

# ------------------------------------------------------------
# Simulation for vector mean problem in R^100 using subsampling
#
# Data:
#   For i=1,...,n_sample:
#     x_i ~ t_5
#     X_i in R^{100} with entries ~ chi^2_1
#     D_i = x_i * X_i   (row-wise scaling)
#
# Parameter:
#   theta_0 = 0 in R^{100}
#
# Estimator:
#   theta_hat = sample mean vector of the observed D_i's
#
# Subsampling:
#   theta_star^(b) = sample mean vector of the subsampled D_i's
#
# Statistics:
#   S_n   = max_{1 <= j <= 100} | sqrt(n) e_j^T (theta_hat - theta_0) |
#         = sqrt(n) * max_j |theta_hat_j|
#
#   S_n^b = max_{1 <= j <= 100} | sqrt(k) e_j^T (theta_star^(b) - theta_hat) |
#
# Usual subsampling coverage check:
#   theta_0 is covered iff S_n <= S_( ceiling(B(1-alpha)) )
#
# Modified subsampling coverage check:
#   theta_0 is covered iff S_n <= S_( ceiling((B+1)(1-alpha)) )
# ------------------------------------------------------------
simulate_subsampling_vector_mean <- function(
    B_grid   = c(15, 20, 25, 50, 100, 150, 200),
    n_iter   = 1000,
    n_sample = 100,
    alpha    = 0.10,
    seed     = 123,
    k        = NULL
) {
  set.seed(seed)

  p <- 100
  theta0 <- rep(0, p)

  if (is.null(k)) {
    k <- max(2L, floor(n_sample^(2/3)))
  } else {
    k <- as.integer(k)
  }

  if (k >= n_sample) {
    stop("k must be strictly smaller than n_sample for subsampling.")
  }

  out <- data.frame(
    B = integer(),
    method = character(),
    coverage = numeric(),
    mc_se = numeric(),
    stringsAsFactors = FALSE
  )

  for (B in B_grid) {
    covered_usual <- logical(n_iter)
    covered_mod   <- logical(n_iter)

    for (iter in seq_len(n_iter)) {
      # ----- Generate D_i as n_sample x p matrix -----
      x_t  <- rt(n_sample, df = 5)
      Xmat <- matrix(rchisq(n_sample * p, df = 1), nrow = n_sample, ncol = p)
      Dmat <- Xmat * x_t

      # ----- theta_hat: sample mean vector in R^p -----
      theta_hat <- colMeans(Dmat)

      # ----- S_n = sqrt(n) * max_j |theta_hat_j - theta0_j| -----
      S_n <- sqrt(n_sample) * max(abs(theta_hat - theta0))

      # ----- Subsampled mean vectors theta_star^(b) -----
      sub_means <- replicate(B, {
        idx <- sample.int(n_sample, size = k, replace = FALSE)
        colMeans(Dmat[idx, , drop = FALSE])
      })

      # Ensure sub_means is p x B even when B = 1
      if (B == 1L) {
        sub_means <- matrix(sub_means, nrow = p, ncol = 1)
      }

      # ----- Subsampling statistics S_n^b -----
      centered_sub <- sub_means - theta_hat
      S_sub <- sqrt(k) * apply(abs(centered_sub), 2, max)
      S_sorted <- sort(S_sub)

      # ----- Usual subsampling -----
      idx_usual <- ceiling(B * (1 - alpha))
      idx_usual <- max(1L, min(B, idx_usual))
      covered_usual[iter] <- (S_n <= S_sorted[idx_usual])

      # ----- Modified subsampling -----
      idx_mod <- ceiling((B + 1) * (1 - alpha))
      idx_mod <- max(1L, min(B, idx_mod))
      covered_mod[iter] <- (S_n <= S_sorted[idx_mod])
    }

    # ----- Summaries -----
    cov_usual <- mean(covered_usual)
    cov_mod   <- mean(covered_mod)

    se_usual <- sqrt(cov_usual * (1 - cov_usual) / n_iter)
    se_mod   <- sqrt(cov_mod   * (1 - cov_mod)   / n_iter)

    out <- rbind(
      out,
      data.frame(B = B, method = "Usual subsampling",    coverage = cov_usual, mc_se = se_usual),
      data.frame(B = B, method = "Modified subsampling", coverage = cov_mod,   mc_se = se_mod)
    )
  }

  out
}

# ------------------------------------------------------------
# Coverage plot only
# ------------------------------------------------------------
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
      title = "Coverage",
      y = "Estimated Coverage",
      x = "Number of Subsamples (B)",
      color = "Method"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 14),
      legend.text  = element_text(size = 13),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      plot.title   = element_text(face = "bold", size = 18, hjust = 0.5)
    )
}

# ------------------------------------------------------------
# Example usage
# ------------------------------------------------------------
res_df_multi <- simulate_subsampling_vector_mean(
  B_grid   = seq(1, 100),
  n_iter   = 1000,
  n_sample = 5000,
  alpha    = 0.10,
  seed     = 123
)

p_cov <- plot_coverage(res_df_multi, alpha = 0.10, n_ticks = 10)
print(p_cov)
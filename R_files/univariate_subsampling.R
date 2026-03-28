############################
###############################

suppressWarnings({
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
})
library(ggplot2)

# ------------------------------------------------------------
# Simulation (parallelized, no plotting)
# ------------------------------------------------------------
simulate_subsampling <- function(
    B_grid   = c(15, 20, 25, 50, 100),
    n_iter   = 1000,
    n_sample = 50,
    alpha    = 0.10,
    seed     = 123,
    k        = NULL,
    n_cores  = parallel::detectCores(logical = TRUE)
) {
  set.seed(seed)

  theta_true <- 1

  if (is.null(k)) {
    k <- max(2L, floor(n_sample^(2 / 3)))
  } else {
    k <- as.integer(k)
  }

  if (k >= n_sample) {
    stop("k must be strictly smaller than n_sample for subsampling.")
  }

  n_cores <- as.integer(n_cores)
  if (is.na(n_cores) || n_cores < 1L) {
    n_cores <- 1L
  }

  # tau_n = n, tau_k = k
  tau_ratio <- k / n_sample

  out <- data.frame(
    B = integer(),
    method = character(),
    coverage = numeric(),
    mc_se = numeric(),
    mean_width = numeric(),
    width_se = numeric(),
    stringsAsFactors = FALSE
  )

  for (B in B_grid) {
    iter_results <- parallel::mclapply(
      X = seq_len(n_iter),
      FUN = function(iter) {
        x <- runif(n_sample, min = 0, max = 1)
        theta_hat <- max(x)

        sub_stats <- replicate(B, {
          idx <- sample.int(n_sample, size = k, replace = FALSE)
          max(x[idx])
        })
        sub_stats_sorted <- sort(sub_stats)

        # ---- Usual subsampling CI ----
        lower_idx_usual <- ceiling(B * (1 - alpha / 2))
        upper_idx_usual <- ceiling(B * (alpha / 2))
        lower_idx_usual <- max(1L, min(B, lower_idx_usual))
        upper_idx_usual <- max(1L, min(B, upper_idx_usual))

        L_usual <- theta_hat + tau_ratio * (theta_hat - sub_stats_sorted[lower_idx_usual])
        U_usual <- theta_hat + tau_ratio * (theta_hat - sub_stats_sorted[upper_idx_usual])

        # ---- Modified subsampling CI ----
        lower_idx_mod <- ceiling((1 - alpha) * (B + 1)) + floor((alpha / 2) * (B + 1))
        upper_idx_mod <- floor((alpha / 2) * (B + 1))
        lower_idx_mod <- max(1L, min(B, lower_idx_mod))
        upper_idx_mod <- max(1L, min(B, upper_idx_mod))

        L_mod <- theta_hat + tau_ratio * (theta_hat - sub_stats_sorted[lower_idx_mod])
        U_mod <- theta_hat + tau_ratio * (theta_hat - sub_stats_sorted[upper_idx_mod])

        list(
          cov_usual = as.integer(theta_true >= L_usual & theta_true < U_usual),
          cov_mod = as.integer(theta_true >= L_mod & theta_true < U_mod),
          width_usual = U_usual - L_usual,
          width_mod = U_mod - L_mod
        )
      },
      mc.cores = n_cores
    )

    covered_usual <- sum(vapply(iter_results, function(x) x$cov_usual, integer(1)))
    covered_mod <- sum(vapply(iter_results, function(x) x$cov_mod, integer(1)))
    widths_usual <- vapply(iter_results, function(x) x$width_usual, numeric(1))
    widths_mod <- vapply(iter_results, function(x) x$width_mod, numeric(1))

    cov_usual <- covered_usual / n_iter
    cov_mod <- covered_mod / n_iter

    se_usual <- sqrt(cov_usual * (1 - cov_usual) / n_iter)
    se_mod <- sqrt(cov_mod * (1 - cov_mod) / n_iter)

    mw_usual <- mean(widths_usual)
    mw_mod <- mean(widths_mod)
    wse_usual <- sd(widths_usual) / sqrt(n_iter)
    wse_mod <- sd(widths_mod) / sqrt(n_iter)

    out <- rbind(
      out,
      data.frame(
        B = B,
        method = "Usual subsampling",
        coverage = cov_usual,
        mc_se = se_usual,
        mean_width = mw_usual,
        width_se = wse_usual
      ),
      data.frame(
        B = B,
        method = "Modified subsampling",
        coverage = cov_mod,
        mc_se = se_mod,
        mean_width = mw_mod,
        width_se = wse_mod
      )
    )
  }

  return(out)
}

# ------------------------------------------------------------
# Plotting helpers
# ------------------------------------------------------------
.make_ticks <- function(B_vals, n_ticks = 10) {
  unique(as.integer(seq(min(B_vals), max(B_vals), length.out = n_ticks)))
}

plot_coverage <- function(df, alpha, n_ticks = 10) {
  tick_breaks <- .make_ticks(df$B, n_ticks)
  ggplot(df, aes(x = B, y = coverage, color = method, group = method)) +
    geom_hline(yintercept = 1 - alpha, linetype = "dotted",
               color = "black", linewidth = 1.5) +
    geom_line() +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = tick_breaks) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(
      title = "Mean Coverage",
      y = "Estimated Coverage",
      x = "Number of Subsamples (B)",
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
      title = "Mean CI Width",
      y = "Mean Width",
      x = "Number of Subsamples (B)",
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

# ------------------------------------------------------------
# Example usage
# ------------------------------------------------------------
res_df_multi <- simulate_subsampling(
  B_grid = seq(1, 200),
  n_iter = 1000,
  n_sample = 100,
  alpha = 0.10,
  seed = 123
)

save(res_df_multi, file = "uniform_max_subsampling.RData")

p_cov <- plot_coverage(res_df_multi, alpha = 0.10, n_ticks = 10)
p_wid <- plot_width(res_df_multi, n_ticks = 10)

dir.create("PDF_images", recursive = TRUE, showWarnings = FALSE)
ggsave(filename = file.path("PDF_images", "univ_sub_cov.pdf"), plot = p_cov, device = "pdf", width = 6, height = 4)
ggsave(filename = file.path("PDF_images", "univ_sub_width.pdf"), plot = p_wid, device = "pdf", width = 6, height = 4)

print(p_cov)
print(p_wid)

suppressWarnings({
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
})
library(ggplot2)

# ------------------------------------------------------------
# Multivariate bootstrap miscoverage simulation in R^100
# - Data dimension p = 100
# - n_sample = 100
# - For each alpha in alpha_grid and c in {2, 4, 6},
#   B(alpha, c) = ceiling(c / alpha - 1)
# - We evaluate realized miscoverage of the modified bootstrap:
#   miscoverage = P(S_n > S_(ceiling((B+1)(1-alpha))))
# ------------------------------------------------------------
simulate_multivar_miscoverage <- function(
    alpha_grid = c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),
    c_grid = c(2, 4, 6),
    n_iter = 1000,
    n_sample = 100,
    p = 100,
    seed = 123,
    n_cores = parallel::detectCores(logical = TRUE)
) {
  set.seed(seed)

  theta0 <- rep(0, p)

  n_cores <- as.integer(n_cores)
  if (is.na(n_cores) || n_cores < 1L) {
    n_cores <- 1L
  }

  out <- data.frame(
    alpha = numeric(),
    B = integer(),
    c_value = integer(),
    curve = character(),
    miscoverage = numeric(),
    mc_se = numeric(),
    ci90_low = numeric(),
    ci90_high = numeric(),
    stringsAsFactors = FALSE
  )

  z90 <- qnorm(0.95)

  for (alpha in alpha_grid) {
    for (c_val in c_grid) {
      B <- ceiling(c_val / alpha - 1)

      covered <- parallel::mclapply(
        X = seq_len(n_iter),
        FUN = function(iter) {
          x_t <- rt(n_sample, df = 5)
          Xmat <- matrix(rchisq(n_sample * p, df = 1), nrow = n_sample, ncol = p)
          Dmat <- Xmat * x_t

          theta_hat <- colMeans(Dmat)
          S_n <- sqrt(n_sample) * max(abs(theta_hat - theta0))

          boot_means <- replicate(B, {
            idx <- sample.int(n_sample, size = n_sample, replace = TRUE)
            colMeans(Dmat[idx, , drop = FALSE])
          })

          if (B == 1L) {
            boot_means <- matrix(boot_means, nrow = p, ncol = 1)
          }

          centered_boot <- boot_means - theta_hat
          S_boot <- sqrt(n_sample) * apply(abs(centered_boot), 2, max)
          S_sorted <- sort(S_boot)

          idx_mod <- ceiling((B + 1) * (1 - alpha))
          idx_mod <- max(1L, min(B, idx_mod))

          as.integer(S_n <= S_sorted[idx_mod])
        },
        mc.cores = n_cores
      )

      covered_vec <- as.integer(unlist(covered))
      cov_hat <- mean(covered_vec)
      mis_hat <- 1 - cov_hat
      se_hat <- sqrt(mis_hat * (1 - mis_hat) / n_iter)
      ci_low <- max(0, mis_hat - z90 * se_hat)
      ci_high <- min(1, mis_hat + z90 * se_hat)

      out <- rbind(
        out,
        data.frame(
          alpha = alpha,
          B = B,
          c_value = c_val,
          curve = paste0("B= (", c_val, "/alpha) - 1"),
          miscoverage = mis_hat,
          mc_se = se_hat,
          ci90_low = ci_low,
          ci90_high = ci_high
        )
      )
    }
  }

  out
}

plot_multivar_miscoverage <- function(df) {
  df$curve <- factor(
    df$curve,
    levels = c(
      "B= (2/alpha) - 1",
      "B= (4/alpha) - 1",
      "B= (6/alpha) - 1"
    )
  )

  ggplot(df, aes(x = alpha, y = miscoverage, color = curve, fill = curve, group = curve)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "red", linewidth = 1.0) +
    geom_ribbon(aes(ymin = ci90_low, ymax = ci90_high), alpha = 0.18, linewidth = 0, show.legend = FALSE) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = sort(unique(df$alpha))) +
    labs(
      title = "Realized Miscoverage vs Nominal Miscoverage",
      x = expression(alpha),
      y = "Realized Miscoverage",
      color = "Bootstrap Size Rule",
      fill = "Bootstrap Size Rule"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.title.y = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )
}

# ------------------------------------------------------------
# Run once, save simulation output, then plot from saved output.
# This allows plot formatting updates without rerunning simulation.
# ------------------------------------------------------------
run_and_save_plot <- function(alpha_grid, seed, rds_file, rdata_file, pdf_file, n_iter = 2500) {
  miscoverage_df <- simulate_multivar_miscoverage(
    alpha_grid = alpha_grid,
    c_grid = c(2, 4, 6),
    n_iter = n_iter,
    n_sample = 100,
    p = 100,
    seed = seed
  )

  saveRDS(miscoverage_df, file = rds_file)
  save(miscoverage_df, file = rdata_file)

  # Plot from stored results so formatting can be changed later
  miscoverage_plot_df <- readRDS(rds_file)
  p_mis <- plot_multivar_miscoverage(miscoverage_plot_df)

  dir.create("PDF_images", recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = file.path("PDF_images", pdf_file),
    plot = p_mis,
    device = "pdf",
    width = 8,
    height = 5
  )

  print(p_mis)
}

# Original alpha sequence
run_and_save_plot(
  alpha_grid = c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),
  seed = 123,
  rds_file = "multivariate_bootstrap_alpha_miscoverage.rds",
  rdata_file = "multivariate_bootstrap_alpha_miscoverage.RData",
  pdf_file = "multivariate_bootstrap_alpha_miscoverage.pdf"
)

# New alpha sequence: 0.05 to 0.95 by 0.05
run_and_save_plot(
  alpha_grid = seq(0.05, 0.95, by = 0.05),
  seed = 456,
  rds_file = "multivariate_bootstrap_alpha_miscoverage_05_to_95.rds",
  rdata_file = "multivariate_bootstrap_alpha_miscoverage_05_to_95.RData",
  pdf_file = "multivariate_bootstrap_alpha_miscoverage_05_to_95.pdf"
)

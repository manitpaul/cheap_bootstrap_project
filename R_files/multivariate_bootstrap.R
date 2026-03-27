suppressWarnings({
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
})
library(ggplot2)

# ------------------------------------------------------------
# Simulation (no plotting)
#   For i=1..n_sample:
#     x_i ~ t_5
#     X_i in R^{100} with entries ~ chi^2_1
#     D_i = x_i * X_i  (elementwise across the 100 coordinates)
#   Statistic:
#     theta_hat = max( colMeans(D_1,...,D_n) )
#   Bootstrap:
#     For each b, resample rows {D_i} with replacement,
#     compute max( colMeans(resample) ).
#   theta_true = 0; coverage checks if 0 in CI (half-open [L,U)).
# ------------------------------------------------------------
simulate_bootstrap <- function(
    B_grid   = c(15, 20, 25, 50, 100),
    n_iter   = 1000,
    n_sample = 50,
    alpha    = 0.10,
    rate     = 5,   # unused; kept for interface compatibility
    seed     = 123
) {
  set.seed(seed)
  theta_true <- 0
  p <- 100   # dimension
  
  out <- data.frame(
    B = integer(), method = character(),
    coverage = numeric(), mc_se = numeric(),
    mean_width = numeric(), width_se = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (B in B_grid) {
    covered_A <- 0L; covered_Bm <- 0L; covered_C <- 0L
    widths_A <- widths_Bm <- widths_C <- numeric(n_iter)
    
    for (iter in seq_len(n_iter)) {
      # ----- Generate D_i as n_sample x p matrix -----
      x_t  <- rt(n_sample, df = 5)                            # t_5 scalars
      Xmat <- matrix(rchisq(n_sample * p, df = 1),            # chi^2_1 entries
                     nrow = n_sample, ncol = p)
      Dmat <- Xmat * x_t                                      # row i scaled by x_i
      
      # ----- Statistic: max coordinate of mean vector -----
      theta_hat <- max(colMeans(Dmat))
      
      # ----- Bootstrap statistics: reapply same statistic to resamples -----
      boot_stats <- replicate(B, {
        idx <- sample.int(n_sample, size = n_sample, replace = TRUE)
        max(colMeans(Dmat[idx, , drop = FALSE]))
      })
      boot_stats_sorted <- sort(boot_stats)
      
      # ---- Usual bootstrap (A): basic / reverse-percentile, ceil/ceil on B ----
      l_idx_A <- ceiling((alpha / 2) * B)
      r_idx_A <- ceiling((1 - alpha / 2) * B)
      l_idx_A <- max(1, min(B, l_idx_A)); r_idx_A <- max(1, min(B, r_idx_A))
      l0_A <- boot_stats_sorted[l_idx_A]; r0_A <- boot_stats_sorted[r_idx_A]
      L_A <- 2 * theta_hat - r0_A; U_A <- 2 * theta_hat - l0_A
      widths_A[iter] <- U_A - L_A
      covered_A <- covered_A + as.integer(theta_true >= L_A & theta_true < U_A)
      
      # ---- Modified bootstrap (B): basic, floor/ceil on (B+1) ----
      l_idx_Bm <- floor((alpha / 2) * (B + 1))
      r_idx_Bm <- ceiling((1 - alpha / 2) * (B + 1))
      l_idx_Bm <- max(1, min(B, l_idx_Bm)); r_idx_Bm <- max(1, min(B, r_idx_Bm))
      l0_Bm <- boot_stats_sorted[l_idx_Bm]; r0_Bm <- boot_stats_sorted[r_idx_Bm]
      L_Bm <- 2 * theta_hat - r0_Bm; U_Bm <- 2 * theta_hat - l0_Bm
      widths_Bm[iter] <- U_Bm - L_Bm
      covered_Bm <- covered_Bm + as.integer(theta_true >= L_Bm & theta_true < U_Bm)
      
      # ---- Cheap bootstrap (C): theta_hat ± t_{B,1-α/2} * S ----
      # S^2 = (1/B) * Σ (theta_b* - theta_hat)^2
      S2 <- mean((boot_stats - theta_hat)^2); S <- sqrt(S2)
      t_q <- qt(1 - alpha/2, df = B)  # 1 - α/2 quantile with B df
      L_C <- theta_hat - t_q * S
      U_C <- theta_hat + t_q * S
      widths_C[iter] <- U_C - L_C
      covered_C <- covered_C + as.integer(theta_true >= L_C & theta_true < U_C)
    }
    
    # Summaries per B
    cov_A  <- covered_A  / n_iter
    cov_Bm <- covered_Bm / n_iter
    cov_C  <- covered_C  / n_iter
    
    se_A  <- sqrt(cov_A  * (1 - cov_A )  / n_iter)
    se_Bm <- sqrt(cov_Bm * (1 - cov_Bm) / n_iter)
    se_C  <- sqrt(cov_C  * (1 - cov_C )  / n_iter)
    
    mw_A  <- mean(widths_A);  wse_A  <- sd(widths_A)  / sqrt(n_iter)
    mw_Bm <- mean(widths_Bm); wse_Bm <- sd(widths_Bm) / sqrt(n_iter)
    mw_C  <- mean(widths_C);  wse_C  <- sd(widths_C)  / sqrt(n_iter)
    
    out <- rbind(
      out,
      data.frame(B=B, method="Usual bootstrap",    coverage=cov_A,  mc_se=se_A,  mean_width=mw_A,  width_se=wse_A),
      data.frame(B=B, method="Modified bootstrap", coverage=cov_Bm, mc_se=se_Bm, mean_width=mw_Bm, width_se=wse_Bm),
      data.frame(B=B, method="Cheap bootstrap",    coverage=cov_C,  mc_se=se_C,  mean_width=mw_C,  width_se=wse_C)
    )
  }
  
  return(out)
}

# ------------------------------------------------------------
# Plotting helpers (outside the simulator)
# ------------------------------------------------------------
.make_ticks <- function(B_vals, n_ticks = 10) {
  unique(as.integer(seq(min(B_vals), max(B_vals), length.out = n_ticks)))
}

plot_coverage <- function(df, alpha, n_ticks = 10) {
  tick_breaks <- .make_ticks(df$B, n_ticks)
  ggplot(df, aes(x = B, y = coverage, color = method, group = method)) +
    geom_hline(yintercept = 1 - alpha, linetype = "dotted", color = "black", linewidth = 1.5) +
    geom_line() + geom_point(size = 1.8) +
    scale_x_continuous(breaks = tick_breaks) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(
      title = "Mean Coverage",
      y = "Estimated Coverage",
      x = "Number of Bootstraps (B)",
      color = "Method"
    ) + theme_bw()
}

plot_width <- function(df, n_ticks = 10) {
  tick_breaks <- .make_ticks(df$B, n_ticks)
  ggplot(df, aes(x = B, y = mean_width, color = method, group = method)) +
    geom_line() + geom_point(size = 1.8) +
    scale_x_continuous(breaks = tick_breaks) +
    labs(
      title = "Mean CI Width",
      y = "Mean Width",
      x = "Number of Bootstraps (B)",
      color = "Method"
    ) + theme_bw()
}

plot_coverage <- function(df, alpha, n_ticks = 10) {
  tick_breaks <- .make_ticks(df$B, n_ticks)
  ggplot(df, aes(x = B, y = coverage, color = method, group = method)) +
    geom_hline(yintercept = 1 - alpha, linetype = "dotted", color = "black", linewidth = 1.5) +
    geom_line() + geom_point(size = 1.8) +
    scale_x_continuous(breaks = tick_breaks) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(
      title = "Mean Coverage",
      y = "Estimated Coverage",
      x = "Number of Bootstraps (B)",
      color = "Method"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 14),
      legend.text  = element_text(size = 13),       # bigger legend labels
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      plot.title   = element_text(face = "bold", size = 18, hjust = 0.5)  # centered & bigger
    )
}

plot_width <- function(df, n_ticks = 10) {
  tick_breaks <- .make_ticks(df$B, n_ticks)
  ggplot(df, aes(x = B, y = mean_width, color = method, group = method)) +
    geom_line() + geom_point(size = 1.8) +
    scale_x_continuous(breaks = tick_breaks) +
    labs(
      title = "Mean CI Width",
      y = "Mean Width",
      x = "Number of Bootstraps (B)",
      color = "Method"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 14),
      legend.text  = element_text(size = 13),       # bigger legend labels
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      plot.title   = element_text(face = "bold", size = 18, hjust = 0.5)  # centered & bigger
    )
}


# ------------------------------------------------------------
# Example usage
# ------------------------------------------------------------
res_df_multi <- simulate_bootstrap(
  B_grid   = seq(1, 200),
  n_iter   = 1000,
  n_sample = 100,
  alpha    = 0.10,
  seed     = 123
)

save(res_df_multi, file = "vector_mean_max_bootstrap.RData")

p_cov <- plot_coverage(res_df_multi, alpha = 0.10, n_ticks = 10)
p_wid <- plot_width(res_df_multi, n_ticks = 10)

print(p_cov)
print(p_wid)

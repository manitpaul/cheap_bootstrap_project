# Percentile bootstrap coverage comparison (two index rules)
# - Data: Exponential(rate = 5), true mean = 1/5
# - Methods:
#   (A) l = ceil((alpha/2) * B),        r = ceil((1 - alpha/2) * B)
#   (B) l = floor((alpha/2) * (B + 1)), r = ceil((1 - alpha/2) * (B + 1))
# - CI form (both): [ 2*theta_hat - r,  2*theta_hat - l )   (basic bootstrap / reverse-percentile)
# - Coverage check uses [lower, upper) as specified

suppressWarnings({
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
})
library(ggplot2)

simulate_coverage <- function(
    B_grid   = c(15, 20, 25, 50, 100),
    n_iter   = 1000,
    n_sample = 50,
    alpha    = 0.10,
    rate     = 5,
    seed     = 123
) {
  set.seed(seed)
  theta_true <- 1 / rate
  
  # storage
  out <- data.frame(
    B = integer(),
    method = character(),
    coverage = numeric(),
    mc_se = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (B in B_grid) {
    covered_A <- 0L
    covered_B <- 0L
    
    for (iter in seq_len(n_iter)) {
      # 1) generate data and estimate mean
      x <- rexp(n_sample, rate = rate)
      theta_hat <- mean(x)
      
      # 2) bootstrap means
      boot_means <- replicate(B, mean(sample(x, size = n_sample, replace = TRUE)))
      boot_means <- sort(boot_means)
      
      # ---------- Method A: ceil, ceil on B ----------
      l_idx_A <- ceiling((alpha / 2) * B)
      r_idx_A <- ceiling((1 - alpha / 2) * B)
      # ensure 1..B (ceil already >=1, but clamp anyway)
      l_idx_A <- max(1, min(B, l_idx_A))
      r_idx_A <- max(1, min(B, r_idx_A))
      
      l0_A <- boot_means[l_idx_A]
      r0_A <- boot_means[r_idx_A]
      
      L_A <- 2 * theta_hat - r0_A
      U_A <- 2 * theta_hat - l0_A
      # [L, U) containment as requested
      covered_A <- covered_A + as.integer(theta_true >= L_A & theta_true < U_A)
      
      # ---------- Method B: floor, ceil on (B+1), clamped to 1..B ----------
      l_idx_B <- floor((alpha / 2) * (B + 1))
      r_idx_B <- ceiling((1 - alpha / 2) * (B + 1))
      l_idx_B <- max(1, min(B, l_idx_B))
      r_idx_B <- max(1, min(B, r_idx_B))
      
      l0_B <- boot_means[l_idx_B]
      r0_B <- boot_means[r_idx_B]
      
      L_B <- 2 * theta_hat - r0_B
      U_B <- 2 * theta_hat - l0_B
      covered_B <- covered_B + as.integer(theta_true >= L_B & theta_true < U_B)
    }
    
    cov_A <- covered_A / n_iter
    cov_B <- covered_B / n_iter
    se_A  <- sqrt(cov_A * (1 - cov_A) / n_iter)
    se_B  <- sqrt(cov_B * (1 - cov_B) / n_iter)
    
    out <- rbind(
      out,
      data.frame(B = B, method = "Rule A: ceil/ceil on B",          coverage = cov_A, mc_se = se_A),
      data.frame(B = B, method = "Rule B: floor/ceil on (B+1)",     coverage = cov_B, mc_se = se_B)
    )
  }
  
  # plot
  p <- ggplot(out, aes(x = B, y = coverage, color = method, group = method)) +
    geom_hline(yintercept = 1 - alpha, linetype = "dotted", color = "red") +
    geom_line() +
    geom_point(size = 2) +
    scale_x_continuous(breaks = B_grid) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = sprintf("Coverage of Percentile Bootstrap CIs (alpha = %.2f, n_iter = %d, n = %d, B varies)", alpha, n_iter, n_sample),
      y = "Estimated Coverage",
      x = "Number of Bootstraps (B)",
      color = "Index Rule"
    ) +
    theme_bw()
  
  list(summary = out, plot = p)
}

# ---- Run with your requested settings ----
res <- simulate_coverage(
  B_grid = seq(1,200),
  n_iter = 1000,
  n_sample = 1000, # change if desired
  alpha = 0.10,
  rate = 5,
  seed = 123
)

res$summary
print(res$plot)

suppressWarnings({
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
})
library(ggplot2)

# Setup 1: Exponential (Rate = 5) Mean estimation
# ------------------------------------------------------------
# Simulation (no plotting)
# ------------------------------------------------------------
simulate_bootstrap <- function(
    B_grid   = c(15, 20, 25, 50, 100),
    n_iter   = 1000,
    n_sample = 50,
    alpha    = 0.10,
    rate     = 5,
    seed     = 123
) {
  set.seed(seed)
  theta_true <- 1 / rate
  
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
      # data & point estimate
      x <- rexp(n_sample, rate = rate)
      theta_hat <- mean(x)
      
      # bootstrap means
      boot_means <- replicate(B, mean(sample(x, size = n_sample, replace = TRUE)))
      boot_means_sorted <- sort(boot_means)
      
      # ---- Usual bootstrap (A) ----
      l_idx_A <- ceiling((alpha / 2) * B)
      r_idx_A <- ceiling((1 - alpha / 2) * B)
      l_idx_A <- max(1, min(B, l_idx_A)); r_idx_A <- max(1, min(B, r_idx_A))
      l0_A <- boot_means_sorted[l_idx_A]; r0_A <- boot_means_sorted[r_idx_A]
      L_A <- 2 * theta_hat - r0_A; U_A <- 2 * theta_hat - l0_A
      widths_A[iter] <- U_A - L_A
      covered_A <- covered_A + as.integer(theta_true >= L_A & theta_true < U_A)
      
      # ---- Modified bootstrap (B) ----
      l_idx_Bm <- floor((alpha / 2) * (B + 1))
      r_idx_Bm <- ceiling((1 - alpha / 2) * (B + 1))
      l_idx_Bm <- max(1, min(B, l_idx_Bm)); r_idx_Bm <- max(1, min(B, r_idx_Bm))
      l0_Bm <- boot_means_sorted[l_idx_Bm]; r0_Bm <- boot_means_sorted[r_idx_Bm]
      L_Bm <- 2 * theta_hat - r0_Bm; U_Bm <- 2 * theta_hat - l0_Bm
      widths_Bm[iter] <- U_Bm - L_Bm
      covered_Bm <- covered_Bm + as.integer(theta_true >= L_Bm & theta_true < U_Bm)
      
      # ---- Cheap bootstrap (C) ----
      # S^2 = (1/B) * sum (theta_b* - theta_hat)^2
      S2 <- mean((boot_means - theta_hat)^2); S <- sqrt(S2)
      t_q <- qt(1 - alpha/2, df = B)  # t_{B, 1-α/2}
      L_C <- theta_hat - t_q * S; U_C <- theta_hat + t_q * S
      widths_C[iter] <- U_C - L_C
      covered_C <- covered_C + as.integer(theta_true >= L_C & theta_true < U_C)
    }
    
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
  
  # return only data (plotting is separate)
  return(out)
}

# ------------------------------------------------------------
# Plotting helpers (outside the simulate function)
# ------------------------------------------------------------
.make_ticks <- function(B_vals, n_ticks = 10) {
  unique(as.integer(seq(min(B_vals), max(B_vals), length.out = n_ticks)))
}

# plot_coverage <- function(df, alpha, n_ticks = 10) {
#   tick_breaks <- .make_ticks(df$B, n_ticks)
#   ggplot(df, aes(x = B, y = coverage, color = method, group = method)) +
#     geom_hline(yintercept = 1 - alpha, linetype = "dotted", color = "black", linewidth = 1.5) +
#     geom_line() + geom_point(size = 1.8) +
#     scale_x_continuous(breaks = tick_breaks) +
#     coord_cartesian(ylim = c(0.5, 1)) +
#     labs(
#       title = "Mean Coverage",
#       y = "Estimated Coverage",
#       x = "Number of Bootstraps (B)",
#       color = "Method"
#     ) + theme_bw()
# }
# 
# plot_width <- function(df, n_ticks = 10) {
#   tick_breaks <- .make_ticks(df$B, n_ticks)
#   ggplot(df, aes(x = B, y = mean_width, color = method, group = method)) +
#     geom_line() + geom_point(size = 1.8) +
#     scale_x_continuous(breaks = tick_breaks) +
#     labs(
#       title = "Mean CI Width",
#       y = "Mean Width",
#       x = "Number of Bootstraps (B)",
#       color = "Method"
#     ) + theme_bw()
# }

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
res_df <- simulate_bootstrap(
  B_grid   = seq(1, 200),
  n_iter   = 1000,
  n_sample = 1000,
  alpha    = 0.10,
  rate     = 5,
  seed     = 123
)

save(res_df, file = "exp_bootstrap.RData") 
p_cov <- plot_coverage(res_df, alpha = 0.10, n_ticks = 10)
p_wid <- plot_width(res_df, n_ticks = 10)

print(p_cov)
print(p_wid)
#############################
############################
###############################

suppressWarnings({
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
})
library(ggplot2)

# ------------------------------------------------------------
# Simulation (no plotting)
#   For each iteration:
#     X_1, ..., X_n ~ U(0,1)
#   Statistic:
#     theta_hat = max(X_1, ..., X_n)
#   Subsampling:
#     For each b, draw a subsample of size k = floor(n^(2/3))
#     without replacement, and compute
#     theta_hat^*_b = max(subsample).
#   Parameter:
#     theta_true = 1
#   Rate:
#     tau_n = n, tau_k = k
# ------------------------------------------------------------
simulate_subsampling <- function(
    B_grid   = c(15, 20, 25, 50, 100),
    n_iter   = 1000,
    n_sample = 50,
    alpha    = 0.10,
    seed     = 123,
    k        = NULL
) {
  set.seed(seed)
  
  theta_true <- 1
  
  if (is.null(k)) {
    k <- max(2L, floor(n_sample^(2/3)))
  } else {
    k <- as.integer(k)
  }
  
  if (k >= n_sample) {
    stop("k must be strictly smaller than n_sample for subsampling.")
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
    covered_usual <- 0L
    covered_mod   <- 0L
    widths_usual  <- numeric(n_iter)
    widths_mod    <- numeric(n_iter)
    
    for (iter in seq_len(n_iter)) {
      # ----- Generate data from U(0,1) -----
      x <- runif(n_sample, min = 0, max = 1)
      
      # ----- Full-sample statistic -----
      theta_hat <- max(x)
      
      # ----- Subsampling statistics -----
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
      
      widths_usual[iter] <- U_usual - L_usual
      covered_usual <- covered_usual +
        as.integer(theta_true >= L_usual & theta_true < U_usual)
      
      # ---- Modified subsampling CI ----
      lower_idx_mod <- ceiling((B + 1) * (1 - alpha / 2))
      upper_idx_mod <- floor((B + 1) * (alpha / 2))
      
      lower_idx_mod <- max(1L, min(B, lower_idx_mod))
      upper_idx_mod <- max(1L, min(B, upper_idx_mod))
      
      L_mod <- theta_hat + tau_ratio * (theta_hat - sub_stats_sorted[lower_idx_mod])
      U_mod <- theta_hat + tau_ratio * (theta_hat - sub_stats_sorted[upper_idx_mod])
      
      widths_mod[iter] <- U_mod - L_mod
      covered_mod <- covered_mod +
        as.integer(theta_true >= L_mod & theta_true < U_mod)
    }
    
    # ----- Summaries per B -----
    cov_usual <- covered_usual / n_iter
    cov_mod   <- covered_mod   / n_iter
    
    se_usual <- sqrt(cov_usual * (1 - cov_usual) / n_iter)
    se_mod   <- sqrt(cov_mod   * (1 - cov_mod)   / n_iter)
    
    mw_usual  <- mean(widths_usual)
    mw_mod    <- mean(widths_mod)
    wse_usual <- sd(widths_usual) / sqrt(n_iter)
    wse_mod   <- sd(widths_mod)   / sqrt(n_iter)
    
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
      legend.text  = element_text(size = 13),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      plot.title   = element_text(face = "bold", size = 18, hjust = 0.5)
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
      legend.text  = element_text(size = 13),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      plot.title   = element_text(face = "bold", size = 18, hjust = 0.5)
    )
}

# ------------------------------------------------------------
# Example usage
# ------------------------------------------------------------
res_df_multi <- simulate_subsampling(
  B_grid   = seq(1, 200),
  n_iter   = 1000,
  n_sample = 100,
  alpha    = 0.10,
  seed     = 123
)

save(res_df_multi, file = "uniform_max_subsampling.RData")

p_cov <- plot_coverage(res_df_multi, alpha = 0.10, n_ticks = 10)
p_wid <- plot_width(res_df_multi, n_ticks = 10)

print(p_cov)
print(p_wid)

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
  
  left_idx  <- ceiling(B * alpha / 2)
  right_idx <- ceiling(B * (1 - alpha / 2))
  
  left_idx  <- max(1L, min(B, left_idx))
  right_idx <- max(1L, min(B, right_idx))
  
  L <- 2 * theta_bar_1 - z[right_idx]
  U <- 2 * theta_bar_1 - z[left_idx]
  
  c(L = L, U = U)
}

modified_sgd_ci <- function(theta_bar_1, theta_star_bar_1, alpha = 0.10) {
  B <- length(theta_star_bar_1)
  z <- sort(theta_star_bar_1)
  
  left_idx  <- floor((B + 1) * alpha / 2)
  right_idx <- ceiling((B + 1) * (1 - alpha / 2))
  
  left_idx  <- max(1L, min(B, left_idx))
  right_idx <- max(1L, min(B, right_idx))
  
  L <- 2 * theta_bar_1 - z[right_idx]
  U <- 2 * theta_bar_1 - z[left_idx]
  
  c(L = L, U = U)
}

# ============================================================
# Main Monte Carlo comparison
# NOTE:
# - The data generator below matches the paper's LAD setting
#   when tau = 0.5.
# - For general tau != 0.5, replace the error generator so that
#   the tau-quantile of epsilon is 0.
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
    seed = 123
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
  
  B_max <- max(B_grid)
  theta01_true <- theta0[1]
  
  covered_usual <- integer(length(B_grid))
  covered_modified <- integer(length(B_grid))
  
  widths_usual <- matrix(NA_real_, nrow = n_iter, ncol = length(B_grid))
  widths_modified <- matrix(NA_real_, nrow = n_iter, ncol = length(B_grid))
  
  for (iter in seq_len(n_iter)) {
    # p = 3
    X <- matrix(rnorm(N * 3), nrow = N, ncol = 3)
    
    # LAD setting: epsilon ~ DE(0,1)
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
      
      widths_usual[iter, j] <- ci_usual["U"] - ci_usual["L"]
      widths_modified[iter, j] <- ci_modified["U"] - ci_modified["L"]
      
      covered_usual[j] <- covered_usual[j] +
        as.integer(theta01_true >= ci_usual["L"] &&
                     theta01_true <  ci_usual["U"])
      
      covered_modified[j] <- covered_modified[j] +
        as.integer(theta01_true >= ci_modified["L"] &&
                     theta01_true <  ci_modified["U"])
    }
  }
  
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
      title = expression(paste("Coverage for CI of ", theta[0,1])),
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
      title = expression(paste("Mean CI width for ", theta[0,1])),
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

print(p_cov)
print(p_wid)

head(res_df)
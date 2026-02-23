compute_SR_one_sample <- function(exp_v, adj_m, maxSR) {
  
  n <- length(exp_v)
  
  # Mass action: sumexp_i = sum_j A_ij * x_j
  sumexp_v <- as.vector(adj_m %*% matrix(exp_v, ncol = 1))
  
  # Stationary distribution: pi_i = x_i * sumexp_i / norm
  invP_v <- exp_v * sumexp_v
  nf <- sum(invP_v)
  
  if (nf == 0) return(list(SR = NA, local_S = rep(NA, n)))
  
  invP_v <- invP_v / nf
  
  # Stochastic matrix: p_ij = A_ij * x_j / sumexp_i
  safe_sumexp <- sumexp_v
  safe_sumexp[safe_sumexp == 0] <- 1
  p_m <- t(t(adj_m) * exp_v) / safe_sumexp
  
  # Local entropy: S_i = -sum_j p_ij * log(p_ij)
  # Vectorized using rowSums
  log_p <- matrix(0, nrow = nrow(p_m), ncol = ncol(p_m))
  nonzero <- p_m > 0
  log_p[nonzero] <- log(p_m[nonzero])
  S_v <- -rowSums(p_m * log_p)
  
  # Entropy rate: SR = sum(pi_i * S_i) / maxSR
  SR <- sum(invP_v * S_v) / maxSR
  
  return(list(SR = SR, local_S = S_v))
}


compute_SR_all_samples <- function(expr_sub, adj_m, maxSR) {
  
  n_samples <- ncol(expr_sub)
  SR_vec <- numeric(n_samples)
  
  for (s in 1:n_samples) {
    exp_v <- as.numeric(expr_sub[, s])
    result <- compute_SR_one_sample(exp_v, adj_m, maxSR)
    SR_vec[s] <- result$SR
  }
  
  return(SR_vec)
}


run_permutation_corrected <- function(adj_m, expr_sub, maxSR,
                                      original_SR,
                                      n_permutations = 100,
                                      seed = 456) {
  set.seed(seed)
  
  n_samples <- ncol(expr_sub)
  n_genes <- nrow(expr_sub)
  
  permuted_SR <- matrix(NA, nrow = n_permutations, ncol = n_samples)
  colnames(permuted_SR) <- colnames(expr_sub)
  
  cat("Running", n_permutations, "permutations (corrected SR)...\n")
  
  for (perm in 1:n_permutations) {
    if (perm %% 10 == 0) cat(sprintf("  Permutation %d/%d\r", perm, n_permutations))
    
    # Shuffle expression within each sample column
    expr_permuted <- expr_sub
    for (s in 1:n_samples) {
      expr_permuted[, s] <- sample(expr_sub[, s])
    }
    
    # Compute SR for all samples with permuted data
    permuted_SR[perm, ] <- compute_SR_all_samples(expr_permuted, adj_m, maxSR)
  }
  
  cat("\n\n")
  
  # Statistics
  p_values <- numeric(n_samples)
  z_scores <- numeric(n_samples)
  
  for (s in 1:n_samples) {
    perm_dist <- permuted_SR[, s]
    obs_val <- original_SR[s]
    
    # Two-tailed p-value
    p_values[s] <- mean(abs(perm_dist - mean(perm_dist)) >= abs(obs_val - mean(perm_dist)))
    
    # Z-score
    z_scores[s] <- (obs_val - mean(perm_dist)) / sd(perm_dist)
  }
  
  cat("=== PERMUTATION RESULTS (Corrected SR) ===\n")
  cat(sprintf("  Samples with p < 0.05: %d / %d (%.1f%%)\n",
              sum(p_values < 0.05), n_samples,
              100 * sum(p_values < 0.05) / n_samples))
  cat(sprintf("  Mean |Z-score|: %.3f\n", mean(abs(z_scores))))
  cat(sprintf("  Z-score range: %.3f to %.3f\n", min(z_scores), max(z_scores)))
  cat(sprintf("  Observed SR mean: %.4f\n", mean(original_SR)))
  cat(sprintf("  Permuted SR mean: %.4f\n\n", mean(colMeans(permuted_SR, na.rm = TRUE))))
  
  if (mean(p_values < 0.05) > 0.8) {
    cat("  -> HIGHLY SIGNIFICANT: Network captures real biological signal\n\n")
  } else if (mean(p_values < 0.05) > 0.5) {
    cat("  -> MODERATELY SIGNIFICANT\n\n")
  } else {
    cat("  -> NOT SIGNIFICANT\n\n")
  }
  
  results <- list(
    permuted_means = permuted_SR,
    original_means = original_SR,
    p_values = p_values,
    z_scores = z_scores
  )
  
  return(results)
}

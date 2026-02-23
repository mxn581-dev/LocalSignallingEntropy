run_perturbation_corrected <- function(network, adj_m, expr_sub, maxSR,
                                       original_SR,
                                       perturbation_rates = c(0.05, 0.10, 0.20),
                                       n_iterations = 100,
                                       seed = 789) {
  set.seed(seed)
  
  n_samples <- ncol(expr_sub)
  mc_genes <- rownames(adj_m)
  original_edges <- ecount(network)
  all_nodes <- V(network)$name
  
  results_list <- list()
  
  for (rate in perturbation_rates) {
    
    cat(sprintf("\n=== Testing %.0f%% edge perturbation (corrected SR) ===\n", rate * 100))
    
    perturbed_SR <- matrix(NA, nrow = n_iterations, ncol = n_samples)
    colnames(perturbed_SR) <- colnames(expr_sub)
    edge_counts <- numeric(n_iterations)
    
    for (iter in 1:n_iterations) {
      if (iter %% 10 == 0) cat(sprintf("  Iteration %d/%d\r", iter, n_iterations))
      
      # ---- Perturb the network ----
      net_pert <- network
      
      n_to_modify <- round(original_edges * rate)
      n_to_delete <- round(n_to_modify / 2)
      n_to_add <- round(n_to_modify / 2)
      
      # Delete random edges
      if (n_to_delete > 0 && ecount(net_pert) > 0) {
        edges_to_delete <- sample(1:ecount(net_pert),
                                  min(n_to_delete, ecount(net_pert)))
        net_pert <- delete_edges(net_pert, edges_to_delete)
      }
      
      # Add random edges (batch)
      if (n_to_add > 0) {
        new_from <- sample(all_nodes, n_to_add * 2, replace = TRUE)
        new_to   <- sample(all_nodes, n_to_add * 2, replace = TRUE)
        keep <- new_from != new_to
        new_from <- new_from[keep]
        new_to   <- new_to[keep]
        if (length(new_from) > n_to_add) {
          new_from <- new_from[1:n_to_add]
          new_to   <- new_to[1:n_to_add]
        }
        if (length(new_from) > 0) {
          edge_vec <- as.vector(rbind(new_from, new_to))
          net_pert <- add_edges(net_pert, edge_vec)
        }
      }
      
      net_pert <- simplify(net_pert)
      edge_counts[iter] <- ecount(net_pert)
      
      # ---- Build perturbed adjacency matrix for mc_genes ----
      # Subset to the same genes used in the original analysis
      genes_in_pert <- intersect(mc_genes, V(net_pert)$name)
      net_sub <- induced_subgraph(net_pert, genes_in_pert)
      
      # Build adjacency matrix aligned to mc_genes
      adj_pert <- matrix(0, nrow = length(mc_genes), ncol = length(mc_genes))
      rownames(adj_pert) <- mc_genes
      colnames(adj_pert) <- mc_genes
      
      # Fill in edges that exist in perturbed network
      adj_sub <- as_adjacency_matrix(net_sub, sparse = FALSE)
      shared <- intersect(mc_genes, rownames(adj_sub))
      adj_pert[shared, shared] <- adj_sub[shared, shared]
      diag(adj_pert) <- 0
      
      # ---- Compute maxSR for perturbed network ----
      # Use the same formula: log(largest eigenvalue)
      fa <- function(x, extra = NULL) { as.vector(adj_pert %*% x) }
      ap <- tryCatch({
        arpack(fa, options = list(n = nrow(adj_pert), nev = 1, which = "LM"), sym = TRUE)
      }, error = function(e) { NULL })
      
      if (is.null(ap) || ap$values[1] <= 1) {
        # Fallback: skip this iteration
        perturbed_SR[iter, ] <- NA
        next
      }
      
      maxSR_pert <- log(ap$values[1])
      
      # ---- Compute corrected SR ----
      perturbed_SR[iter, ] <- compute_SR_all_samples(expr_sub, adj_pert, maxSR_pert)
    }
    
    cat("\n")
    
    # Statistics
    cv_per_sample <- apply(perturbed_SR, 2, sd, na.rm = TRUE) /
                     colMeans(perturbed_SR, na.rm = TRUE)
    overall_cv <- mean(cv_per_sample, na.rm = TRUE)
    mean_shift <- mean(colMeans(perturbed_SR, na.rm = TRUE) - original_SR, na.rm = TRUE)
    
    cat(sprintf("  Overall CV: %.4f (%.2f%%)\n", overall_cv, overall_cv * 100))
    cat(sprintf("  Mean edges: %.0f (original: %d, change: %+.1f%%)\n",
                mean(edge_counts), original_edges,
                100 * (mean(edge_counts) - original_edges) / original_edges))
    cat(sprintf("  Mean SR shift from original: %+.4f\n", mean_shift))
    
    if (overall_cv < 0.01) {
      cat("  -> VERY ROBUST\n")
    } else if (overall_cv < 0.05) {
      cat("  -> ROBUST\n")
    } else if (overall_cv < 0.10) {
      cat("  -> MODERATE\n")
    } else {
      cat("  -> SENSITIVE\n")
    }
    
    results_list[[paste0("rate_", rate)]] <- list(
      perturbation_rate = rate,
      perturbed_means = perturbed_SR,
      original_means = original_SR,
      cv_per_sample = cv_per_sample,
      overall_cv = overall_cv,
      mean_shift = mean_shift,
      edge_counts = edge_counts
    )
  }
  
  cat("\n===== PERTURBATION SUMMARY =====\n")
  for (rate_name in names(results_list)) {
    r <- results_list[[rate_name]]
    cat(sprintf("  %s: CV = %.4f (%.2f%%), shift = %+.4f\n",
                rate_name, r$overall_cv, r$overall_cv * 100, r$mean_shift))
  }
  cat("================================\n\n")
  
  return(results_list)
}
#   adj_m         - adjacency matrix
#   maxSR         - maximum entropy rate
#   mc_genes      - genes in maximally connected component
#   expr_corrected - linear-scale expression matrix
#   SR_per_sample  - original SR values

# Subset expression to mc_genes
expr_sub <- expr_corrected[mc_genes, , drop = FALSE]

perm_results <- run_permutation_corrected(
  adj_m = adj_m,
  expr_sub = expr_sub,
  maxSR = maxSR,
  original_SR = SR_per_sample,
  n_permutations = 100,
  seed = 456
)

# Export
perm_summary <- data.frame(
  sample = colnames(expr_sub),
  original_SR = perm_results$original_means,
  permuted_mean = colMeans(perm_results$permuted_means, na.rm = TRUE),
  permuted_sd = apply(perm_results$permuted_means, 2, sd, na.rm = TRUE),
  z_score = perm_results$z_scores,
  p_value = perm_results$p_values
)
write.csv(perm_summary, "permutation_corrected_summary.csv", row.names = FALSE)
write.csv(perm_results$permuted_means, "permutation_corrected_means.csv", row.names = FALSE)
cat("Saved: permutation_corrected_summary.csv, permutation_corrected_means.csv\n\n")


perturb_results <- run_perturbation_corrected(
  network = ppi_net,
  adj_m = adj_m,
  expr_sub = expr_sub,
  maxSR = maxSR,
  original_SR = SR_per_sample,
  perturbation_rates = c(0.05, 0.10, 0.20),
  n_iterations = 100,
  seed = 789
)

# Export
for (rate_name in names(perturb_results)) {
  result <- perturb_results[[rate_name]]
  
  write.csv(result$perturbed_means,
            paste0("perturbation_corrected_", rate_name, "_means.csv"),
            row.names = FALSE)
  
  write.csv(data.frame(
    sample = colnames(expr_sub),
    cv = result$cv_per_sample,
    original_SR = result$original_means,
    perturbed_mean = colMeans(result$perturbed_means, na.rm = TRUE)
  ), paste0("perturbation_corrected_", rate_name, "_summary.csv"), row.names = FALSE)

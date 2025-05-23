# Function to calculate residuals from linear regression
calculate_residuals <- function(X, i, U_i_j) {
  # If U_i_j is empty, return the original variable
  if (length(U_i_j) == 0) {
    return(scale(X[, i], center = TRUE, scale = FALSE))
  }
  
  # Create design matrix for regression
  X_U <- X[, U_i_j, drop = FALSE]
  
  # Fit linear regression
  model <- lm(X[, i] ~ X_U)
  
  # Return residuals
  return(residuals(model))
}

# Function to calculate the normalized partial correlation test statistic
calculate_partial_correlation <- function(X, i, j, U_i_j) {
  # First try direct linear regression
  e_i <- calculate_residuals(X, i, U_i_j)
  e_j <- calculate_residuals(X, j, U_i_j)
  
  # Convert to numeric vectors immediately
  e_i <- as.numeric(e_i)
  e_j <- as.numeric(e_j)
  
  # Check if residuals have matching dimensions
  if (length(e_i) != length(e_j)) {
    # Use lasso regression with glmnet for dimension mismatch case
    require(glmnet)
    
    # Center and scale the response variables
    y_i <- scale(X[, i], center = TRUE, scale = FALSE)
    y_j <- scale(X[, j], center = TRUE, scale = FALSE)
    
    # Center and scale the shared neighbors
    X_U <- scale(X[, U_i_j, drop = FALSE], center = TRUE, scale = TRUE)
    
    # Use cross-validation to find optimal lambda for i
    cv_fit_i <- cv.glmnet(X_U, y_i, alpha = 1, nfolds = 5)
    fit_i <- glmnet(X_U, y_i, alpha = 1, lambda = cv_fit_i$lambda.min)
    e_i <- as.numeric(y_i - predict(fit_i, X_U))
    
    # Use cross-validation to find optimal lambda for j
    cv_fit_j <- cv.glmnet(X_U, y_j, alpha = 1, nfolds = 5)
    fit_j <- glmnet(X_U, y_j, alpha = 1, lambda = cv_fit_j$lambda.min)
    e_j <- as.numeric(y_j - predict(fit_j, X_U))
  }
  
  # Calculate sample covariance of residuals
  hat_r_ij <- cov(e_i, e_j)

  # Calculate sample covariance of X_i and X_j
  hat_S_ij_X <- cov(X[, i], X[, j])
  # Calculate D_ij
  D_ij <- hat_S_ij_X - hat_r_ij
  
  # Calculate empirical variance estimator
  hat_theta_ij <- mean((e_i * e_j - hat_r_ij)^2)
  
  # Calculate normalized partial correlation test statistic
  hat_T_ij <- sqrt(nrow(X)) * hat_r_ij / sqrt(hat_theta_ij)
  
  return(list(hat_T_ij = hat_T_ij, hat_theta_ij = hat_theta_ij, D_ij = D_ij, hat_r_ij = hat_r_ij))
}


find_union_neighbors <- function(i, j, non_zero_idx) {
  # First, ensure we're working with the full edge set, not just upper triangular
  # Create a complete edge set by adding reversed edges
  full_edges <- rbind(non_zero_idx, non_zero_idx[, c(2, 1), drop = FALSE])
  
  # Get all neighbors of i (nodes connected to i)
  edges_with_i <- full_edges[full_edges[, 1] == i, , drop = FALSE]
  neighbors_i <- unique(edges_with_i[, 2])
  
  # Get all neighbors of j (nodes connected to j)
  edges_with_j <- full_edges[full_edges[, 1] == j, , drop = FALSE]
  neighbors_j <- unique(edges_with_j[, 2])
  
  # Find union of neighbors excluding i and j
  U_i_j <- setdiff(union(neighbors_i, neighbors_j), c(i, j))
  
  return(U_i_j)
}


# Main function to get selected pairs using partial correlation
get_selected_pairs_partial <- function(alpha, 
                                      data, 
                                      non_zero_idx_upper, 
                                      verbose = TRUE) {
  p <- ncol(data)
  
  # Calculate Bonferroni-corrected probability
  prob_value <- alpha / (2 * p^2)
  
  # Ensure probability is strictly within (0, 1) for qnorm
  prob_value <- min(max(prob_value, 1e-10), 1 - 1e-10)
  
  # Calculate threshold
  gamma_ij <- qnorm(1 - prob_value)
  
  # Initialize selected pairs
  selected_idx <- c()
  
  # For each edge in the initial estimated graph
  for (edge_idx in 1:nrow(non_zero_idx_upper)) {
    i <- non_zero_idx_upper[edge_idx, 1]
    j <- non_zero_idx_upper[edge_idx, 2]
    
    # Find union of neighbors
    U_i_j <- find_union_neighbors(i, j, non_zero_idx_upper)
    
    # Calculate partial correlation test statistic
    pc_result <- calculate_partial_correlation(data, i, j, U_i_j)
    hat_T_ij <- pc_result$hat_T_ij
    
    # Check if edge should be included based on partial correlation
    if (abs(hat_T_ij) >= gamma_ij) {
      selected_idx <- c(selected_idx, edge_idx)
    }
  }
  
  # Return selected pairs
  return(non_zero_idx_upper[selected_idx, , drop = FALSE])
}



# Function to get FDR and TPR remains the same
get_FDR_TPR <- function(selected_pairs, adj_mat_True){
    # Get dimensions from true adjacency matrix
    p <- nrow(adj_mat_True)
    
    # Create estimated adjacency matrix
    adj_mat_est <- matrix(0, p, p)
    
    # Fill in the estimated adjacency matrix
    if (nrow(selected_pairs) > 0) {
      for (i in 1:nrow(selected_pairs)) {
        adj_mat_est[selected_pairs[i, 1], selected_pairs[i, 2]] <- 1
      }
    }
    
    # Make matrix symmetric
    adj_mat_est <- adj_mat_est + t(adj_mat_est)
    
    # Set diagonals to 1
    diag(adj_mat_est) <- 1
    
    # Convert to regular matrix if needed
    adj_mat_True <- as.matrix(adj_mat_True)
    
    # Calculate true positives and false positives
    true_neg <- which(adj_mat_True == 0, arr.ind = TRUE)
    true_pos <- which(adj_mat_True == 1, arr.ind = TRUE)
    
    # False positives: should be 0 but marked as 1
    FP <- sum(adj_mat_est[true_neg])
    
    # True positives: should be 1 and marked as 1 (excluding diagonal)
    TP <- sum(adj_mat_est[true_pos]) - p
    
    # Calculate total positive edges in true matrix (excluding diagonal)
    pos_sum <- sum(adj_mat_True[true_pos]) - p
    
    # Calculate TPR and FDR
    TPR <- TP / pos_sum
    FDR <- if(FP + TP > 0) FP / (FP + TP) else NA

    adj_mat_est <- as.numeric(adj_mat_est)
    adj_mat_True <- as.numeric(adj_mat_True)
    
    # get hamming distance  
    hamming_distance <- sum(abs(adj_mat_True - adj_mat_est)) / 2
    
    # get MSE by glasso with zero pattern
    S <- cov(data)
    zeroes <- which(adj_mat_est == 0, arr.ind = TRUE)
    est_glasso <- glasso(S, rho = 0, zero = zeroes, nobs = n,
      thr = tol / 100,
      approx = FALSE,
      trace = 0,
      penalize.diagonal = FALSE
    )
    loglik_glasso <- -est_glasso$loglik
    
    return(list(TPR = TPR, FDR = FDR, hamming_distance = hamming_distance, loglik_glasso = loglik_glasso))
}

# Example usage - test the algorithm with different alpha values
test_partial_correlation_algorithm <- function(data, non_zero_idx_upper, adj_mat_True, alpha_range, method, savefile=TRUE) {
  results <- lapply(alpha_range, function(alpha) {
    # Select pairs using partial correlation
    print(alpha)
    selected_pairs <- get_selected_pairs_partial(alpha, data, non_zero_idx_upper)
    
    # Calculate FDR and TPR
    metrics <- get_FDR_TPR(selected_pairs, adj_mat_True)
    
    df <- data.frame(
      method = method,
      idx = NA,
      time = NA,
      alpha = alpha,
      loglik = metrics$loglik_glasso,
      hd_cri = metrics$hamming_distance,
      mse_cri = -1,
      tpr_percentage = metrics$TPR,
      fdr_percentage = metrics$FDR
    )
    return(df)
  })
  
  # Convert to data frame
  results_df <- do.call(rbind, lapply(results, function(x) as.data.frame(as.list(x))))
  
  if (savefile){
    save_path <- paste0(temp_save_dir, "partial_corr_results_",method, "_", alpha_range[1], "_", alpha_range[length(alpha_range)], "_", gtype, "_", model, "_", n, "_", p, "_", nt, ".csv")
    cat("CSV saved to: ", save_path, "\n")
    write.csv(results_df, file = save_path, row.names = FALSE)
  }
  return(results_df)
}
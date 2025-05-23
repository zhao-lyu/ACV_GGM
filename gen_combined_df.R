#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(glasso)
source(paste0(getwd(), "/metric.R"))

# Function to read a specific result file
read_result_file <- function(g, m, n, p, nt, method, base_dir) {
  # Construct the filename
  filename <- paste0("results_", g, "_", m, "_", n, "_", p, "_", nt, "_", method, ".rds")
  filepath <- file.path(base_dir, filename)
  
  # Check if file exists
  if (!file.exists(filepath)) {
    warning(paste("File not found:", filepath))
    return(NULL)
  }
  
  # Read the RDS file
  results <- readRDS(filepath)
  return(results)
}

# Function to calculate metrics for a result
calculate_metrics <- function(result, data, adj_mat_True, theta_0, fixed_lam, tol = 1e-3) {
  if (is.null(result)) return(NULL)
  
  # Get the optimal index
  opt_idx <- result$idx
  
  # Get the estimation at the optimal index
  if (result$method == "BIC" || result$method == "EBIC") {
    # For BIC, we have the refitted estimates
    est <- result$info$est_w_refitting[[opt_idx]]
  } else {
    # For other methods, we need to compute the estimate at the optimal lambda
    S <- cov(data)
    
    # Debug information
    cat("Method:", result$method, "\n")
    cat("Optimal index:", opt_idx, "\n")
    cat("Length of fixed_lam:", length(fixed_lam), "\n")
    
    # Ensure we're using a single lambda value
    if (opt_idx > length(fixed_lam)) {
      cat("Warning: opt_idx exceeds length of fixed_lam. Using last lambda value.\n")
      lambda_value <- fixed_lam[length(fixed_lam)]
    } else {
      lambda_value <- fixed_lam[opt_idx]
    }
    
    # Check if lambda_value is a scalar
    if (length(lambda_value) > 1) {
      cat("Warning: lambda_value is a vector. Using first element.\n")
      lambda_value <- lambda_value[1]
    }
    
    # Ensure lambda_value is a scalar
    if (length(lambda_value) != 1) {
      cat("Error: lambda_value is not a scalar. Using default value.\n")
      lambda_value <- 0.1  # Default value
    }
    
    cat("Using lambda value:", lambda_value, "\n")
    
    # Call glasso with the scalar lambda value
    est <- glasso(S, rho = lambda_value, 
                 thr = tol/100, approx = FALSE, trace = 0, penalize.diagonal = FALSE)
  }
  
  # Calculate metrics
  hd_cri <- hd_loss(1, list(est), adj_mat_True, use_tol = TRUE, tol = tol)
  mse_cri <- mse_loss(1, list(est), theta_0)
  tpr_percentage <- TPR(1, list(est), adj_mat_True, use_tol = TRUE, tol = tol)
  fdr_percentage <- FDR(1, list(est), adj_mat_True, use_tol = TRUE, tol = tol)
  
  # Calculate negative log-likelihood
  S <- cov(data)
  N <- nrow(data)
  K <- est$wi
  d <- det(K)
  negloglik <- -(N/2) * (log(d) - sum(diag(S %*% K)))
  
  return(list(
    hd_cri = hd_cri,
    mse_cri = mse_cri,
    tpr_percentage = tpr_percentage,
    fdr_percentage = fdr_percentage,
    negloglik = negloglik,
    time = result$total_time
  ))
}

# Function to combine all results into a single data frame
combine_all_results <- function(n_list, p_list, g_list, m_list, nt = 42, base_dir) {
  # Initialize empty lists to store results
  all_results <- list()
  
  # Define methods
  methods <- c("K3CV", "K5CV", "K10CV", "LOOCV", "ACV", "QUIC", "BIC", "EBIC", "ACV_Kfold")

  
  grids <- expand.grid(g = g_list, m = m_list, n = n_list, p = p_list, nt = nt)
  for (grid_idx in 1:nrow(grids)) {
    g <- grids$g[grid_idx]
    m <- grids$m[grid_idx]
    n <- grids$n[grid_idx]
    p <- grids$p[grid_idx]
    nt <- grids$nt[grid_idx]
    
    cur_dir <- getwd()
    src_data_path <- paste0(cur_dir, "/Data/")
    src_g_file <- paste0(src_data_path, g, "_", p, "_", nt, ".rds")
    src_lam_file <- paste0(src_data_path, g, "_", p, "_", nt, "_lamALL.rds")
    
    # Check if data files exist
    if (!file.exists(src_g_file) || !file.exists(src_lam_file)) {
      warning(paste("Data files not found for g=", g, "p=", p, "nt=", nt))
      next
    }
    
    G <- readRDS(src_g_file)
    lambda.all_algo <- readRDS(src_lam_file)
    lambda.all <- lambda.all_algo[[m]]
    
    # Get data and parameters
    adj_mat_True <- G$adj_mat_True
    theta_0 <- G$theta_0
    data <- G$data[1:n,]
    data <- scale(data, center = TRUE, scale = FALSE)
    
    # Get lambda values
    n_list_lambda <- as.integer(rownames(lambda.all))
    n_idx <- which(n_list_lambda == n)
    fixed_lam <- lambda.all[n_idx,]

    print(length(methods))
    
    # Initialize vectors to store metrics for all methods
    tpr_values <- rep(NA, length(methods))
    shd_values <- rep(NA, length(methods))
    fdr_values <- rep(NA, length(methods))
    mse_values <- rep(NA, length(methods))
    negloglik_values <- rep(NA, length(methods))
    time_values <- rep(NA, length(methods))
    
    # check how many methods in methods are processed
    print(paste0("Processing grid index: ", grid_idx))

    # Process each method
    for (i in seq_along(methods)) {
      method <- methods[i]
      results <- read_result_file(g, m, n, p, nt, method, base_dir)
      
      if (!is.null(results)) {
        metrics <- calculate_metrics(results, data, adj_mat_True, theta_0, fixed_lam)
        
        if (!is.null(metrics)) {
          tpr_values[i] <- metrics$tpr_percentage
          shd_values[i] <- metrics$hd_cri
          fdr_values[i] <- metrics$fdr_percentage
          mse_values[i] <- metrics$mse_cri
          negloglik_values[i] <- metrics$negloglik
          time_values[i] <- metrics$time
        }
      }
    }
    
    # Create a data frame for metrics with replicated values
    metrics_df <- data.frame(
      graph = rep(g, length(methods)),
      model = rep(m, length(methods)),
      n = rep(n, length(methods)),
      p = rep(p, length(methods)),
      nt = rep(nt, length(methods)),
      method = methods,
      TPR = tpr_values,
      SHD = shd_values,
      FDR = fdr_values,
      MSE = mse_values,
      negloglik = negloglik_values,
      time = time_values
    )
    
    all_results[[length(all_results) + 1]] <- metrics_df
  }
  
  # Combine all results
  if (length(all_results) > 0) {
    combined_df <- do.call(rbind, all_results)
    return(combined_df)
  } else {
    warning("No results found for the specified parameters.")
    return(NULL)
  }
}

# Main function
main <- function() {
  p_val <- Sys.getenv('P_VAL')
  n_val <- Sys.getenv('N_VAL')
  p_list <- as.integer(strsplit(p_val, ",")[[1]])
  n_list <- as.integer(strsplit(n_val, ",")[[1]])
  g_list <- c("knn","sf","er","band")
  m_list <- c("glasso")
  nt <- c(1)
  base_dir <- "acv_save"
  save_dir <- "df_plots"
  save_name <- paste0("combined_result_p", p_list, ".csv")

  
  # Print the parameters
  cat("Generating combined data frame with the following parameters:\n")
  cat("Sample sizes (n):", paste(n_list, collapse = ", "), "\n")
  cat("Feature sizes (p):", paste(p_list, collapse = ", "), "\n")
  cat("Graph types (g):", paste(g_list, collapse = ", "), "\n")
  cat("Models (m):", paste(m_list, collapse = ", "), "\n")
  cat("Number of trials (nt):", nt, "\n")
  
  # Create save directory if it doesn't exist
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  # Combine all results
  cat("\nCombining results...\n")
  all_results_df <- combine_all_results(n_list, p_list, g_list, m_list, nt, base_dir)
  
  if (is.null(all_results_df)) {
    cat("No results found for the specified parameters.\n")
    return(NULL)
  }
  
  cat("Results combined successfully. Total rows:", nrow(all_results_df), "\n")
  
  # Save the combined data frame to CSV
  output_file <- file.path(save_dir, save_name)
  cat("\nSaving combined data frame to:", output_file, "\n")
  write.csv(all_results_df, output_file, row.names = FALSE)
  
  cat("\nDone!\n")
}

# Run the main function
main() 
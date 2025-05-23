library("R.utils")
library("glasso")
library("flare")

options(mc.allow.recursive = TRUE)
library(parallel)

CV_partition <- function(N, nfold){
  # set seed
  set.seed(123)
  ntest <- floor(N/nfold)
  ntrain <- N - ntest
  # ith col is data used for train/test for ith fold
  train_mat <- matrix(NA, nrow = ntrain, ncol = nfold)
  test_mat <- matrix(NA, nrow = ntest, ncol = nfold)
  # fill in data
  N_seq <- c(1:N)
  data_idx <- sample(N)
  for(i in 1:nfold){
    test_idx <- ((i-1)*ntest+1):(i*ntest)
    train_idx <- N_seq[!N_seq %in% test_idx]
    test_mat[,i] <- data_idx[test_idx]
    train_mat[,i] <- data_idx[train_idx]
  }
  cv_data_list <- list(train_mat=train_mat, test_mat=test_mat)
  return(cv_data_list)
}

# K-fold cross validation for graphical lasso
# Splits data into K folds, trains on K-1 folds and validates on remaining fold
# Repeats for each fold and averages results
# 
# Args:
#   N: Number of samples
#   nfold: Number of folds for cross-validation 
#   data: Input data matrix
#   fixed_lam: Vector of lambda values to test
#   tol: Tolerance threshold for convergence
#
# Returns:
#   List containing:
#   - cv.mean: Mean cross-validation loss for each lambda
#   - cv.sd: Standard deviation of CV loss for each lambda
#   - total_time: Total computation time
CV_KCV <- function(N, nfold, data, fixed_lam, ncores, tol = 1e-3) {
  library(parallel)
  library(glasso)
  
  Nlam <- length(fixed_lam)
  p <- ncol(data)
  cv_data_list <- CV_partition(N, nfold)

  # Precompute training covariance matrices
  S_train_list <- vector("list", nfold)
  for (i in 1:nfold) {
    S_train_list[[i]] <- cov(data[cv_data_list$train_mat[, i], ])
  }

  # Create list of all (fold, lambda) index pairs
  idx_list <- expand.grid(i = 1:nfold, j = 1:Nlam)

  # Parallelized loop over (i, j) pairs
  result_list <- mclapply(1:nrow(idx_list), function(k) {
    i <- idx_list$i[k]
    j <- idx_list$j[k]

    S_train <- S_train_list[[i]]
    S_test <- cov(data[cv_data_list$test_mat[, i], ])

    start_time <- Sys.time()
    res <- glasso(
      s = S_train,
      rho = fixed_lam[j],
      thr = tol / 100,
      approx = FALSE,
      trace = 0,
      penalize.diagonal = FALSE
    )
    end_time <- Sys.time()

    K <- res$wi
    d <- det(K)
    loss <- -(N / 2) * (log(d) - sum(diag(S_test %*% K)))
    time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))

    list(loss = loss, time = time_taken, fold = i, lam_idx = j)
  }, mc.cores = ncores)

  # Reshape result into loss matrix [nfold × Nlam]
  loss_mat <- matrix(NA, nrow = nfold, ncol = Nlam)
  time_mat <- matrix(NA, nrow = nfold, ncol = Nlam)
  for (res in result_list) {
    loss_mat[res$fold, res$lam_idx] <- res$loss
    time_mat[res$fold, res$lam_idx] <- res$time
  }


  cv_loss.mean <- colMeans(loss_mat)
  cv_loss.sd <- apply(loss_mat, 2, sd)
  total_time <- sum(time_mat)
  

  return(list(
    cv.mean = cv_loss.mean,
    cv.sd = cv_loss.sd,
    total_time = total_time
  ))
}

# Leave-One-Out Cross Validation for graphical lasso
# Performs leave-one-out CV where each sample is used as a test set once
# while the remaining samples form the training set
# 
# Args:
#   N: Number of samples
#   data: Input data matrix
#   fixed_lam: Vector of lambda values to test
#   ncores: Number of CPU cores to use for parallel computation
#   tol: Tolerance threshold for convergence (default: 1e-3)
#
# Returns:
#   List containing:
#   - cv.mean: Mean cross-validation loss for each lambda
#   - cv.sd: Standard deviation of CV loss for each lambda
#   - total_time: Total computation time
CV_LOOCV <- function(N, data, fixed_lam, ncores, tol=1e-3) {
  library(parallel)
  Nlam <- length(fixed_lam)
  p <- ncol(data)
  n <- nrow(data)
  S_all <- cov(data)

  # Precompute S_train and S_test for each i
  S_train_list <- vector("list", N)
  S_test_list <- vector("list", N)
  for (i in 1:N) {
    data.test <- data[i, ]
    S_test_list[[i]] <- tcrossprod(data.test)
    S_train_list[[i]] <- cov(data[-i, ])
  }

  total_time <- 0
  start_all <- Sys.time()

  # Create all combinations of i and j
  idx_list <- expand.grid(i = 1:N, j = 1:Nlam)
  # Parallel computation of losses
  results <- mclapply(1:nrow(idx_list), function(k) {
    i <- idx_list$i[k]
    j <- idx_list$j[k]
    
    S_train <- S_train_list[[i]]
    S_test <- S_test_list[[i]]
    start_time <- Sys.time()

    res <- glasso::glasso(
      s = S_train, 
      rho = fixed_lam[j], 
      thr = tol / 100,
      approx = FALSE, 
      trace = 0, 
      penalize.diagonal = FALSE
    )
    end_time <- Sys.time()
    time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))

    K <- res$wi
    d <- det(K)
    loss <- -(N/2) * (log(d) - sum(diag(S_test %*% K)))
    # list(loss = loss, time = time_taken, fold = i, lam_idx = j)
    list(loss = loss, time = time_taken)
  }, mc.cores = ncores)

  loss_mat_flat <- sapply(results, `[[`, "loss")
  time_vec <- sapply(results, `[[`, "time")
  total_time <- sum(time_vec)
  loss_mat <- matrix(loss_mat_flat, nrow = N, ncol = Nlam)
  cv_loss.mean <- colMeans(loss_mat)
  cv_loss.sd <- apply(loss_mat, 2, sd)

  return(list(cv.mean = cv_loss.mean, cv.sd = cv_loss.sd, total_time = total_time))
}

# Approximate Cross Validation for graphical lasso
# Uses first-order approximation to estimate leave-one-out CV loss
# without refitting the model for each fold
# 
# Args:
#   N: Number of samples
#   data: Input data matrix
#   fixed_lam: Vector of lambda values to test
#   ncores: Number of CPU cores to use for parallel computation
#   tol: Tolerance threshold for convergence (default: 1e-3)
#   method: Method to use for approximation (default: "NS")
#
# Returns:
#   List containing:
#   - cv.mean: Mean approximate CV loss for each lambda
#   - cv.sd: Standard deviation of approximate CV loss for each lambda
#   - total_time: Total computation time
CV_ACV <- function(N, data, fixed_lam, ncores, tol = 1e-3, method = "NS") {
  library(glasso)
  library(parallel)
  
  Nlam <- length(fixed_lam)
  S_all <- cov(data)
  n <- nrow(data)
  p <- ncol(data)
  full_est <- glasso::glassopath(s = S_all, rholist = fixed_lam, thr = tol/100, approx = FALSE, trace = 0, penalize.diagonal = FALSE)
  
  # Precompute S_diff and S_test
  S_test_list <- vector("list", N)
  S_diff_list <- vector("list", N)
  for (i in 1:N) {
    data.train <- data[-i, ]
    data.test <- data[i, ]
    S_train <- cov(data.train)
    S_diff_list[[i]] <- S_all - S_train * ((n - 2) / (n - 1))
    S_test_list[[i]] <- tcrossprod(data.test)
  }

  idx_list <- expand.grid(i = 1:N, j = 1:Nlam)
  
  results <- mclapply(1:nrow(idx_list), function(k) {
    i <- idx_list$i[k]
    j <- idx_list$j[k]

    S_diff <- S_diff_list[[i]]
    S_test <- S_test_list[[i]]
    support <- abs(full_est$wi[, , j]) > tol

    S_masked <- matrix(0, nrow = p, ncol = p)
    S_masked[support] <- S_diff[support]

    update <- matrix(0, nrow = p, ncol = p)
    start_time <- Sys.time()
    update_temp <- full_est$wi[, , j] %*% S_masked %*% full_est$wi[, , j]
    update <- matrix(0, nrow = p, ncol = p)
    update[support] <- update_temp[support]

    K <- full_est$wi[, , j] + update
    end_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    loss <- -(N/2) * (log(det(K)) - sum(diag(S_test %*% K)))
    if(is.nan(loss)) {
      loss <- Inf
    }
    list(loss = loss, time = elapsed_time)
  }, mc.cores = ncores)

  loss_mat_flat <- sapply(results, `[[`, "loss")
  time_vec <- sapply(results, `[[`, "time")
  total_time <- sum(time_vec)
  loss_mat <- matrix(loss_mat_flat, nrow = N, ncol = Nlam)
  cv_loss.mean <- colMeans(loss_mat)
  cv_loss.sd <- apply(loss_mat, 2, sd)
  
  return(list(cv.mean = cv_loss.mean, cv.sd = cv_loss.sd, total_time = total_time))
}

# ProxACV for graphical lasso via QUIC
# 
# Args:
#   N: Number of samples
#   data: Input data matrix
#   fixed_lam: Vector of lambda values to test
#   ncores: Number of CPU cores to use for parallel computation
#   tol: Tolerance threshold for convergence (default: 1e-3)
#   active_step: Number of active steps for QUIC algorithm (default: 1)
#
# Returns:
#   List containing:
#   - cv.mean: Mean CV loss for each lambda
#   - cv.sd: Standard deviation of CV loss for each lambda
#   - total_time: Total computation time
CV_QUIC<- function(N, data, fixed_lam, ncores, tol = 1e-3, active_step = 1) {
  library(QUIC)
  library(parallel)
  Nlam <- length(fixed_lam)
  S_all <- cov(data)
  n <- nrow(data)
  p <- ncol(data)

  full_est <- glasso::glassopath(
    s = S_all,
    rholist = fixed_lam,
    thr = tol / 100,
    approx = FALSE,
    trace = 0,
    penalize.diagonal = FALSE
  )

  S_test_list <- vector("list", N)
  S_train_list <- vector("list", N)
  for (i in 1:N) {
    x_i <- data[i, , drop = FALSE]
    S_test_list[[i]] <- t(x_i) %*% x_i
    S_train_list[[i]] <- S_all - S_test_list[[i]] / n
  }

  idx_list <- expand.grid(i = 1:N, j = 1:Nlam)

  results <- mclapply(1:nrow(idx_list), function(k) {
    i <- idx_list$i[k]
    j <- idx_list$j[k]

    S_test <- S_test_list[[i]]
    S_train <- S_train_list[[i]]

    start_time <- Sys.time()
    QUIC_res <- QUIC::QUIC(
      S_train,
      rho = fixed_lam[j],
      msg = 0,
      maxIter = active_step,
      X.init = full_est$wi[, , j],
      W.init = full_est$w[, , j]
    )
    end_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    K <- QUIC_res$X
    detK <- det(K)
    loss <- -(N / 2) * (log(detK) - sum(diag(S_test %*% K)))
    # list(loss = loss, time = elapsed_time, fold = i, lam_idx = j)
    list(loss = loss, time = elapsed_time)
  }, mc.cores = ncores)

  loss_mat_flat <- sapply(results, `[[`, "loss")
  time_vec <- sapply(results, `[[`, "time")
  total_time <- sum(time_vec)
  loss_mat <- matrix(loss_mat_flat, nrow = N, ncol = Nlam)
  cv_loss.mean <- colMeans(loss_mat)
  cv_loss.sd <- apply(loss_mat, 2, sd)

  return(list(cv.mean = cv_loss.mean, cv.sd = cv_loss.sd, total_time = total_time))
}


# BIC/EBIC for graphical lasso
# 
# Args:
#   N: Number of samples
#   data: Input data matrix
#   fixed_lam: Vector of lambda values to test
#   ncores: Number of CPU cores to use for parallel computation
#   tol: Tolerance threshold for convergence (default: 1e-3)
#
# Returns:
#   List containing:
#   - cv.mean: Mean cross-validation loss for each lambda
#   - cv.sd: Standard deviation of CV loss for each lambda
#   - total_time: Total computation time
BIC <- function(N, data, fixed_lam, ncores=ncores, tol = 1e-3, countDiagonal = FALSE, use_tol = TRUE, type = "BIC") {
  library(parallel)
  library(glasso)
  Nlam <- length(fixed_lam)
  S_all <- cov(data)

  # Decide on gamma
  gamma <- if (type == "BIC") 0 else 1

  start_time <- Sys.time()

  results <- mclapply(1:Nlam, function(j) {
    tryCatch({
      rho <- fixed_lam[j]
      
      # Fit with penalty
      res <- glasso(
        s = S_all,
        rho = rho,
        thr = tol / 100,
        approx = FALSE,
        trace = 0,
        penalize.diagonal = FALSE
      )

      # Build adjacency
      if (use_tol) 
      {
        adj_mat_est <- abs(res$wi) > tol
      } 
      else 
      {
        adj_mat_est <- res$wi != 0
      }
      diag(adj_mat_est) <- 1

      # Identify zero pattern
      zeroes <- which(adj_mat_est == 0, arr.ind = TRUE)

      # Refit with fixed sparsity pattern
      res_refit <- glasso(
        S_all,
        rho = 0,
        zero = zeroes,
        nobs = N,
        thr = tol / 100,
        approx = FALSE,
        trace = 0,
        penalize.diagonal = FALSE
      )

      # Count edges
      if (use_tol) {
        E <- sum(abs(res_refit$wi[lower.tri(res_refit$wi, diag = countDiagonal)]) > tol)
      } else {
        E <- sum(res_refit$wi[lower.tri(res_refit$wi, diag = countDiagonal)] != 0)
      }

      # BIC or EBIC loss
      bic_val <- -2 * res_refit$loglik + E * log(N) + 4 * gamma * E * log(nrow(res_refit$wi))
      
      list(bic = bic_val, refit = res_refit)
    }, error = function(e) {
      # Return a list with Inf for BIC value to indicate failure
      list(bic = Inf, refit = NULL)
    })
  }, mc.cores = ncores)

  # Check if any results are NULL and replace with a default value
  results <- lapply(results, function(x) {
    if (is.null(x)) {
      list(bic = Inf, refit = NULL)
    } else {
      x
    }
  })

  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Extract BIC values, replacing NULL with Inf
  bic_values <- sapply(results, function(x) if(is.null(x$bic)) Inf else x$bic)

  return(list(
    loss = bic_values,
    est_w_refitting = lapply(results, `[[`, "refit"),
    total_time = total_time
  ))
}


ACV_Kfold <- function(N, nfold, data, fixed_lam, ncores, tol = 1e-3, method = "NS") {
  library(glasso)
  library(parallel)
  
  Nlam <- length(fixed_lam)
  S_all <- cov(data)
  n <- nrow(data)
  p <- ncol(data)
  
  # Get the CV partition
  cv_data_list <- CV_partition(N, nfold)
  
  # Run glasso on full data for all lambda values
  full_est <- glasso::glassopath(s = S_all, rholist = fixed_lam, thr = tol/100, 
                                approx = FALSE, trace = 0, penalize.diagonal = FALSE)
  
  # Precompute S_diff and S_test for each fold
  S_test_list <- vector("list", nfold)
  S_diff_list <- vector("list", nfold)
  ntest <- floor(N/nfold)
  
  for (i in 1:nfold) {
    test_idx <- cv_data_list$test_mat[, i]
    train_idx <- cv_data_list$train_mat[, i]
    
    data.train <- data[train_idx, ]
    data.test <- data[test_idx, ]
    
    S_train <- cov(data.train)
    S_test_list[[i]] <- cov(data.test)
    
    # Adjust the S_diff calculation for K-fold
    ntrain <- length(train_idx)
    S_diff_list[[i]] <- S_all - S_train * ((ntrain) / (n))
  }
  
  # Create list of all (fold, lambda) index pairs
  idx_list <- expand.grid(i = 1:nfold, j = 1:Nlam)
  
  # Parallelized loop over (i, j) pairs
  results <- mclapply(1:nrow(idx_list), function(k) {
    i <- idx_list$i[k]
    j <- idx_list$j[k]
    
    S_diff <- S_diff_list[[i]]
    S_test <- S_test_list[[i]]
    
    support <- abs(full_est$wi[, , j]) > tol
    S_masked <- matrix(0, nrow = p, ncol = p)
    S_masked[support] <- S_diff[support]
    
    start_time <- Sys.time()
    
    # First-order approximation update
    update_temp <- full_est$wi[, , j] %*% S_masked %*% full_est$wi[, , j]
    update <- matrix(0, nrow = p, ncol = p)
    update[support] <- update_temp[support]
    
    # Approximated precision matrix for this fold and lambda
    K <- full_est$wi[, , j] + update
    
    end_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Calculate loss using test data from this fold
    ntest <- nrow(data[cv_data_list$test_mat[, i], ])
    loss <- -(ntest/2) * (log(det(K)) - sum(diag(S_test %*% K)))
    
    if(is.nan(loss)) {
      loss <- Inf
    }
    
    list(loss = loss, time = elapsed_time, fold = i, lam_idx = j)
  }, mc.cores = ncores)
  
  # Extract results into matrices
  loss_mat <- matrix(NA, nrow = nfold, ncol = Nlam)
  time_mat <- matrix(NA, nrow = nfold, ncol = Nlam)
  
  for (res in results) {
    loss_mat[res$fold, res$lam_idx] <- res$loss
    time_mat[res$fold, res$lam_idx] <- res$time
  }
  
  total_time <- sum(time_mat)
  cv_loss.mean <- colMeans(loss_mat)
  cv_loss.sd <- apply(loss_mat, 2, sd)
  
  return(list(cv.mean = cv_loss.mean, cv.sd = cv_loss.sd, total_time = total_time))
}
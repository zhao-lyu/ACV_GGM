# IMPORT
library("glasso")
source(paste0(getwd(), "/metric.R"))
source(paste0(getwd(), "/fdr_utils.R"))


args <- commandArgs(trailingOnly = TRUE)

# Initialize variables
n <- NULL
p <- NULL
nt <- NULL
model <- NULL
gtype <- NULL

# Parse command line arguments
for (arg in args) {
    split_arg <- strsplit(arg, "=")[[1]]
    switch(split_arg[1],
            "--n" = {n <- as.integer(split_arg[2])},
            "--p" = {p <- as.integer(split_arg[2])},
            "--nt" = {nt <- as.integer(split_arg[2])},
            "--m" = {model <- split_arg[2]},
            "--g" = {gtype <- split_arg[2]})
}

# Check if all required arguments are provided
if (is.null(n) || is.null(p) || is.null(nt) || is.null(model) || is.null(gtype)) {
    warning("Missing required arguments! Please provide all of: --n, --p, --nt, --m, --g")
    stop()
}

# LOAD DATA and File
cur_dir <- getwd()
src_data_path <- paste0(cur_dir, "/Data/")
src_g_file <- paste0(src_data_path,gtype,"_", p, "_", nt, ".rds")
src_lam_file <- paste0(src_data_path,gtype,"_",p,"_", nt, "_lamALL.rds")
cat(src_g_file,"\n")
cat(src_lam_file,"\n")

# GET TEMP SAVE
base_dir <- paste0(cur_dir,"/acv_save")
temp_save_dir <- paste0(cur_dir,"/EXP_FDR/")
if (!dir.exists(temp_save_dir)){
    dir.create(temp_save_dir)
}
temp_name <- paste0("results_", gtype,"_",model,"_",n,"_",p,"_",nt)
output_text <- paste0(temp_save_dir, "test_", gtype,"_",model,"_", n,"_", p,"_",nt, "_FDR.txt")
cat("SAVE TO :", output_text,"\n")

# SET PARAM
use_tol <- TRUE
tol <- 1e-3
refitting <- FALSE

# LOAD DATA and lambda
cat("Loading ", src_g_file, src_lam_file, "............ ")
G <- readRDS(src_g_file)
lambda.all_algo <- readRDS(src_lam_file)
lambda.all <- lambda.all_algo[[model]]
cat("DONE!\n")

# GET DATA
adj_mat_True <- G$adj_mat_True
theta_0 <- G$theta_0
data <- G$data[1:n,]
data <- scale(data, center = TRUE, scale = FALSE)
n_list <- as.integer(rownames(lambda.all))
n_idx <- which(n_list == n)
fixed_lam <- lambda.all[n_idx,]

# HERE IS AN EXAMPLE FOR METHODS TO APPLY POST CLEANING
methods <- c("BIC", "K3CV", "ACV_Kfold")
cat("\nmethods: ", methods, file=output_text, append=TRUE)

all_results <- lapply(methods, function(method){
    filename <- paste0("results_", gtype, "_", model, "_", n, "_", p, "_", nt, "_", method, ".rds")
    filepath <- file.path(base_dir, filename)
    results <- readRDS(filepath)
    return(results)
})

all_idx <- lapply(all_results, function(results){
    return(results$idx)
})
all_idx <- unlist(all_idx)

all_time <- lapply(all_results, function(results){
    return(results$total_time)
})
all_time <- unlist(all_time)

cat("\nall_idx: ", all_idx, file=output_text, append=TRUE)
cat("\nall_time: ", all_time, file=output_text, append=TRUE)


S <- cov(data)
N <- nrow(data)

est_wo_refitting <- lapply(seq_along(all_idx), function(i){
    glasso(S, rho = fixed_lam[all_idx[i]],
        thr = tol/100, approx = FALSE, trace = 0, penalize.diagonal = FALSE)
})


negative_loglik_wo_refitting <- sapply(seq_along(all_idx), function(i){
    K <- est_wo_refitting[[i]]$wi
    # check if K is NULL
    if (is.null(K)){
        cat("ERROR: K is NULL", file=temp_name, append=TRUE)
        return(0)
    }
    d <- det(K)
    -(N/2)*(log(d) - sum(diag(S%*%K)))
})
loglik_wo_refitting <- negative_loglik_wo_refitting
cat("\nloglik_all wo refitting: ", loglik_wo_refitting, file=output_text, append=TRUE)

# with refitting
est_refitting <- lapply(seq_along(all_idx), function(i){
    K <- est_wo_refitting[[i]]$wi
    if (use_tol)
    {
      adj_mat_est <- abs(K) > tol
    }
    else
    {
      adj_mat_est <- K != 0
    }
    diag(adj_mat_est) <- 1
    # Identify zero pattern
    zeroes <- which(adj_mat_est == 0, arr.ind = TRUE)
    # Refit with fixed sparsity pattern
    glasso(
      S,
      rho = 0,
      zero = zeroes,
      nobs = N,
      thr = tol / 100,
      approx = FALSE,
      trace = 0,
      penalize.diagonal = FALSE
    )
})

loglik_refitting <- sapply(seq_along(all_idx), function(i){
    -est_refitting[[i]]$loglik
})



if (refitting){
    est <- est_refitting
    loglik_all <- loglik_refitting
} else {
    est <- est_wo_refitting
    loglik_all <- loglik_wo_refitting
}

hd_cri <- sapply(seq_along(all_idx), function(i){hd_loss(i, est, adj_mat_True)})
mse_cri <-  sapply(seq_along(all_idx), function(i){mse_loss(i, est, theta_0)})
tpr_percentage <- sapply(seq_along(all_idx), function(i){TPR(i, est, adj_mat_True)})
fdr_percentage <- sapply(seq_along(all_idx), function(i){FDR(i, est, adj_mat_True)})
# make a dataframe of all results
no_control_df <- data.frame(
    method = methods,
    idx = all_idx,
    time = all_time,
    alpha = -1,
    loglik = loglik_all,
    hd_cri = hd_cri,
    mse_cri = mse_cri,
    tpr_percentage = tpr_percentage,
    fdr_percentage = fdr_percentage
)

# cat all results to test.txt
cat("\nhd_cri: ", hd_cri, file=output_text, append=TRUE)
cat("\nmse_cri: ", mse_cri, file=output_text, append=TRUE)
cat("\ntpr_percentage: ", tpr_percentage, file=output_text, append=TRUE)
cat("\nfdr_percentage: ", fdr_percentage, file=output_text, append=TRUE)

# define a function to wrap the code below
get_non_zero_idx <- function(pred_mat, use_tol, tol){
    if (use_tol){
        pred_mat <- abs(pred_mat) > tol
    } else {
        pred_mat <- pred_mat != 0
    }

    # Make matrix symmetric using AND operator
    pred_mat <- pred_mat & t(pred_mat)

    # Get non-zero indices
    non_zero_idx <- which(pred_mat != 0, arr.ind = TRUE)

    # Get upper triangular indices
    non_zero_idx_upper <- non_zero_idx[non_zero_idx[,1] < non_zero_idx[,2], , drop = FALSE]
    return(non_zero_idx_upper)
}

library(parallel)

# --- Setup ---
alpha_range <- seq(0, 1.0, 0.1)
param_grid <- expand.grid(method = methods,
                          stringsAsFactors = FALSE)

# --- Parallel-safe wrapper ---
parallel_eval <- function(params) {
  method <- params$method
  
  idx <- which(methods == method)
  one_est <- est[[idx]]
  
  non_zero_idx_upper <- get_non_zero_idx(one_est$wi, use_tol, tol)
  non_zero_idx_upper_true <- get_non_zero_idx(G$adj_mat_True, FALSE, 0)
  
  res <- test_partial_correlation_algorithm(
    data,
    non_zero_idx_upper,
    G$adj_mat_True,
    alpha_range = alpha_range,
    method = method,
    savefile = TRUE
  )
  return(res)
}

# --- Run in parallel ---
results_list <- mclapply(
  split(param_grid, seq_len(nrow(param_grid))),
  parallel_eval,
  mc.cores = detectCores()
)

# --- Combine all results ---
all_res <- do.call(rbind, results_list)

all_results_df <- rbind(no_control_df, all_res)
# sort by method, alpha
all_results_df <- all_results_df[order(all_results_df$method, all_results_df$alpha),]
# save all results to a csv file
write.csv(all_results_df, file = paste0(temp_save_dir, "all_results_", gtype, "_", model, "_", n, "_", p, "_", nt, ".csv"), row.names = FALSE)




cat("SAVE TO :", output_text, file=output_text, append=TRUE)
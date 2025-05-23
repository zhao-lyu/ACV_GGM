library("glasso")
library("flare")
source(paste0(getwd(), "/utils.R"))
source(paste0(getwd(), "/metric.R"))

run_experiment <- function(n, p, nt, model, gtype, method) {
    # Print the arguments
    cat(n,"\n")
    cat(p, "\n")
    cat(nt, "\n")
    cat(model, "\n")
    cat(gtype, "\n")
    cat("Method:", method, "\n")

    # LOAD DATA and File
    cur_dir <- getwd()
    src_data_path <- paste0(cur_dir, "/Data/")
    src_g_file <- paste0(src_data_path,gtype,"_", p, "_", nt, ".rds")
    src_lam_file <- paste0(src_data_path,gtype,"_",p,"_", nt, "_lamALL.rds")
    cat(src_g_file,"\n")
    cat(src_lam_file,"\n")

    # GET TEMP SAVE
    temp_save_dir <- paste0(cur_dir,"/acv_save/")
    temp_txt_save_dir <- paste0(cur_dir,"/acv_save_txt/")
    # create temp_save_dir if not exists
    if (!dir.exists(temp_save_dir)) {
        dir.create(temp_save_dir)
    }
    if (!dir.exists(temp_txt_save_dir)) {
        dir.create(temp_txt_save_dir)
    }
    temp_name <- paste0("results_", gtype,"_",model,"_",n,"_",p,"_",nt,"_",method)
    dst_filename <- paste0(temp_save_dir, temp_name,".rds")
    temp_name_txt <- paste0(temp_txt_save_dir, "test_", gtype,"_",model,"_", p,"_",nt, "_",method,".txt")
    cat("SAVE TO :", dst_filename,"\n")

    # SET PARAM
    use_tol <- TRUE
    tol <- 1e-3
    refitting <- FALSE
    ncores <- 16

    # LOAD DATA and lambda
    cat("Loading ", src_g_file, src_lam_file, "............ ")
    G <- readRDS(src_g_file)
    lambda.all_algo <- readRDS(src_lam_file)
    lambda.all <- lambda.all_algo[[model]]
    cat("DONE!\n")

    # GET DATA 
    adj_mat_True = G$adj_mat_True 
    theta_0 = G$theta_0
    data <- G$data[1:n,]
    data <- scale(data, center = TRUE, scale = FALSE)
    n_list <- as.integer(rownames(lambda.all))
    n_idx <- which(n_list == n)
    fixed_lam <- lambda.all[n_idx,]

    # Run selected methods
    if (method == "K3CV") {
        K3CV_result <- CV_KCV(n, 3, data, fixed_lam, ncores)
        idx_list <- which(K3CV_result$cv.mean == min(K3CV_result$cv.mean))
        opt_idx <- idx_list[length(idx_list)]
        results <- list(
            method = "K3CV",
            info = K3CV_result,
            idx = opt_idx,
            total_time = K3CV_result$total_time
        )
    }
    if (method == "K5CV") {
        K5CV_result <- CV_KCV(n, 5, data, fixed_lam, ncores)
        idx_list <- which(K5CV_result$cv.mean == min(K5CV_result$cv.mean))
        opt_idx <- idx_list[length(idx_list)]
        results <- list(
            method = "K5CV",
            info = K5CV_result,
            idx = opt_idx,
            total_time = K5CV_result$total_time
        )
    }
    if (method == "K10CV") {
        K10CV_result <- CV_KCV(n, 10, data, fixed_lam, ncores)
        idx_list <- which(K10CV_result$cv.mean == min(K10CV_result$cv.mean))
        opt_idx <- idx_list[length(idx_list)]
        results <- list(
            method = "K10CV",
            info = K10CV_result,
            idx = opt_idx,
            total_time = K10CV_result$total_time
        )
    }
    if (method == "LOOCV") {
        LOOCV_result <- CV_LOOCV(n, data, fixed_lam, ncores)
        idx_list <- which(LOOCV_result$cv.mean == min(LOOCV_result$cv.mean))
        opt_idx <- idx_list[length(idx_list)]
        results <- list(
            method = "LOOCV",
            info = LOOCV_result,
            idx = opt_idx,
            total_time = LOOCV_result$total_time
        )
    }
    if (method == "ACV") {
        ACV_result <- CV_ACV(n, data, fixed_lam, ncores)
        idx_list <- which(ACV_result$cv.mean == min(ACV_result$cv.mean))
        opt_idx <- idx_list[length(idx_list)]
        results <- list(
            method = "ACV",
            info = ACV_result,
            idx = opt_idx,
            total_time = ACV_result$total_time
        )
    }
    if (method == "QUIC") {
        Quic_result <- CV_QUIC(n, data, fixed_lam, ncores)
        idx_list <- which(Quic_result$cv.mean == min(Quic_result$cv.mean))
        opt_idx <- idx_list[length(idx_list)]
        results <- list(
            method = "QUIC",
            info = Quic_result,
            idx = opt_idx,
            total_time = Quic_result$total_time
        )
    }
    if (method == "BIC") {
        BIC_result <- BIC(n, data, fixed_lam, ncores)
        idx_list <- which(BIC_result$loss == min(BIC_result$loss))
        opt_idx <- idx_list[length(idx_list)]
        results <- list(
            method = "BIC",
            info = BIC_result,
            idx = opt_idx,
            total_time = BIC_result$total_time
        )
    }

    if (method == "EBIC") {
        EBIC_result <- BIC(n, data, fixed_lam, ncores,type = "EBIC")
        idx_list <- which(EBIC_result$loss == min(EBIC_result$loss))
        opt_idx <- idx_list[length(idx_list)]
        results <- list(
            method = "EBIC",
            info = EBIC_result,
            idx = opt_idx,
            total_time = EBIC_result$total_time
        )
    }

    if (method == "ACV_Kfold") {
        ACV_Kfold_result <- ACV_Kfold(n, 10, data, fixed_lam, ncores)
        idx_list <- which(ACV_Kfold_result$cv.mean == min(ACV_Kfold_result$cv.mean))
        opt_idx <- idx_list[length(idx_list)]
        results <- list(
            method = "ACV_Kfold",
            info = ACV_Kfold_result,
            idx = opt_idx,
            total_time = ACV_Kfold_result$total_time
        )
    }

    cat("\nDone for n = ", n, "p = ", p, "graph = ", gtype, "nt = ", nt, file=temp_name_txt, append=TRUE)
    cat("\nidx: ", results$idx, file=temp_name_txt, append=TRUE)
    cat("\ntime: ", results$total_time, file=temp_name_txt, append=TRUE)


    # save results
    saveRDS(results, dst_filename)
}

# Get command line arguments and run if called directly
if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    
    # Initialize variables
    n <- NULL
    p <- NULL
    nt <- NULL
    model <- NULL
    gtype <- NULL
    method <- "all"  # Default to running all methods
    
    # Parse command line arguments
    for (arg in args) {
        split_arg <- strsplit(arg, "=")[[1]]
        switch(split_arg[1],
                "--n" = {n <- as.integer(split_arg[2])},
                "--p" = {p <- as.integer(split_arg[2])},
                "--nt" = {nt <- as.integer(split_arg[2])},
                "--m" = {model <- split_arg[2]},
                "--g" = {gtype <- split_arg[2]},
                "--method" = {method <- split_arg[2]})
    }
    
    arg_len <- sapply(c(n,p,gtype,model), length)
    na_check <- sapply(c(n,p,gtype,model), is.na)
    if(sum(arg_len) < 4 || any(na_check)){
        warning("Incomplete args!")
        stop()
    }
    
    run_experiment(n, p, nt, model, gtype, method)
}

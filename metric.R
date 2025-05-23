# param: 
#      idx: between 1 to Nlam
#      adj_mat_True: the true adjacency matrix based on theta_0
#      use_tol: take tolerance into consideration to decide the adjaceny mat based on invSigma
# return: the number of mismatch between true adj_mat and adj_mat based on the optimal lasso with lam = fixed_lam[idx]
hd_loss <- function(idx, est, adj_mat_True, use_tol=TRUE, tol=1e-3){
    invSigma <- est[[idx]]$wi
    if(use_tol){
        adj_mat <- ifelse(abs(invSigma) > tol, 1, 0)
    }else{
        adj_mat <- ifelse(invSigma !=0, 1, 0)
    }
    diag(adj_mat) <- 1
    adj_mat_True <- as.matrix(adj_mat_True) # Convert S4 matrix to base R matrix
    sum(abs(adj_mat_True - adj_mat))/2
}


# param: 
#      idx: between 1 to Nlam
#      use_tol: take tolerance into consideration to decide the adjaceny mat based on invSigma
# return: False Positive Rate(FDR)
FDR <- function(idx, est, adj_mat_True, use_tol=TRUE, tol=1e-3,countDiagonal=FALSE, file="err_msg.txt"){
    K <- est[[idx]]$wi
    if(use_tol){
        adj_mat <- ifelse(abs(K) > tol, 1, 0)
    }else{
        adj_mat <- ifelse(K !=0, 1, 0)
    }
    diag(adj_mat) <- 1
    true_neg <- which(adj_mat_True == 0, arr.ind = TRUE)
    true_pos <- which(adj_mat_True == 1, arr.ind = TRUE)
    # should be 0, but marked 1 (FP)
    FP = sum(adj_mat[true_neg])
    # should be 1, and marked 1(TP)
    TP = sum(adj_mat[true_pos])-nrow(adj_mat_True)
    if(FP+TP <= 0){
        cat(paste0("\nAt lambda idx = ", idx, ",  estimated adj mat is all zeros; invalid FDR.\n"), file=file,append=TRUE)
        return(NA)
    }
    return(FP/(FP+TP))
}

# param: 
#      idx: between 1 to Nlam
#      use_tol: take tolerance into consideration to decide the adjaceny mat based on invSigma
# return: Positive Selection Rate (the number of edges in estimated adj_mat/total number of edges of true adj_mat)
TPR <- function(idx, est, adj_mat_True, use_tol=TRUE, tol=1e-3,countDiagonal=FALSE){
    K <- est[[idx]]$wi
    if(use_tol){
        adj_mat <- ifelse(abs(K) > tol, 1, 0)
    }else{
        adj_mat <- ifelse(K !=0, 1, 0)
    }
    diag(adj_mat) <- 1
    true_pos <- which(adj_mat_True == 1, arr.ind = TRUE)
    # should be 1, and marked 1(TP)
    TP <- sum(adj_mat[true_pos])-nrow(adj_mat_True)
    # all pos
    pos_sum <- sum(adj_mat_True[true_pos])-nrow(adj_mat_True)
    # return(sum(adj_mat[true_pos])/total_edge_num)
    return(TP/pos_sum)
}

# param: 
#      idx: between 1 to Nlam
#      theta_0: true precision mat
# return: squared difference between estimated precision mat and true precision mat
mse_loss <- function(idx, est, theta_0){
    invSigma <- est[[idx]]$wi
    sum(theta_0 - invSigma)^2
}

cal_num_edge <- function(K, total_edge_num, tol=1e-3){
    adj_mat_K <- ifelse(abs(K) > tol, 1, 0)
    diag(adj_mat_K) <- 1
    return(sum(adj_mat_K) - nrow(K))
}
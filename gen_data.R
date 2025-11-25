# Generate data based on graph type and P

# generator func
generator <- function(N,P,g,ss,band=1,dist_type="gaussian"){
    set.seed(ss)
    cat("ss: ", ss,"\n")
    # band and sf by flare
    if(g == "band"){
        library('flare')
        G <- sugm.generator(n=N, d=P, graph=g, v=0.3,u=0,g=band,seed=ss)
        theta_0 <- G$omega
        sigma_0 <- G$sigma
        data <- G$data
        S <- G$sigmahat
        adj_mat_True <- G$theta
    }else if(g == "sf"){
        library('flare')
        G <- sugm.generator(n=N, d=P, graph="scale-free",seed=ss)
        theta_0 <- G$omega
        sigma_0 <- G$sigma
        data <- G$data
        S <- G$sigmahat
        adj_mat_True <- G$theta
    }
    else if(g == "er"){
        library("igraph")
        gnp_m <- ceiling(P * log(P) * 0.1)
        G_renyi <- erdos.renyi.game(n = P, p.or.m = gnp_m, type = "gnm", directed = FALSE, loops = FALSE)
        adj_mat_True <- as_adjacency_matrix(graph = G_renyi, type = "both", sparse = TRUE)
        diag(adj_mat_True) <- 0
        v <- 0.3
        u <- 0
        theta_0 <- adj_mat_True*v
        diag(theta_0) <- abs(min(eigen(theta_0)$values)) + 0.1 + u
        sigma_0 <- cov2cor(solve(theta_0))
        theta_0 <- solve(sigma_0)

    }else if(g == "knn"){
        library("mstknnclust")
        library("stats")
        library("igraph")
        # random data
        x <- matrix(runif(P*P, min = -10, max = 10), nrow=P, ncol=P)
        d <- base::as.matrix(stats::dist(x, method="euclidean"))
        cg <- generate.complete.graph(1:nrow(x),d)
        G_knn <- generate.knn(cg, suggested.k=1)
        adj_mat_True <- as_adjacency_matrix(graph = G_knn$knn.graph, type = "both", sparse = TRUE)
        adj_mat_True[which(as.matrix(adj_mat_True) > 0, arr.ind = TRUE)] <- 1
        diag(adj_mat_True) <- 0
        v <- 0.3
        u <- 0
        theta_0 <- adj_mat_True*v
        diag(theta_0) <- abs(min(eigen(theta_0)$values)) + 0.1 + u
        sigma_0 <- cov2cor(solve(theta_0))
        theta_0 <- solve(sigma_0)
    }
    # CONSTRUCT DIST FOR DATA
    library("MASS")
    # mean
    mu <- rep(0, P)
    # Generate data based on specified distribution
    if(dist_type == "gaussian") { # Gaussian
        data <- mvrnorm(N, mu, sigma_0)
    } else if(dist_type == "subt") { # t-distribution
        df <- 30
        Z <- matrix(rt(N * P, df), nrow = N, ncol = P)
        L <- chol(sigma_0)
        data <- Z %*% L
    } else if(dist_type == "subuniform") { # Uniform
        Z <- matrix(runif(N * P, -sqrt(3), sqrt(3)), nrow = N, ncol = P)
        L <- chol(sigma_0)
        data <- Z %*% L
    } else if(dist_type == "subbounded") { # Bounded
        Z <- matrix(2 * rbinom(N * P, 1, 0.5) - 1, nrow = N, ncol = P)
        L <- chol(sigma_0)
        data <- Z %*% L
    } else if(dist_type == "submixture") { # Mixture
        mixture_prob <- 0.7
        mean1 <- rep(0, P)
        mean2 <- rep(3, P)
        
        component <- rbinom(N, 1, mixture_prob)
        
        # Initialize data matrix
        data <- matrix(0, nrow = N, ncol = P)
        
        # Generate from component 1
        n1 <- sum(component == 1)
        if (n1 > 0) {
            Z1 <- matrix(rnorm(n1 * P), nrow = n1, ncol = P)
            L <- chol(sigma_0)
            data[component == 1, ] <- Z1 %*% L + matrix(rep(mean1, each = n1), nrow = n1)
        }
        
        # Generate from component 2
        n2 <- sum(component == 0)
        if (n2 > 0) {
            Z2 <- matrix(rnorm(n2 * P), nrow = n2, ncol = P)
            L <- chol(sigma_0)
            data[component == 0, ] <- Z2 %*% L + matrix(rep(mean2, each = n2), nrow = n2)
        }
    } else if(dist_type == "laplace") { # Laplace
        library(rmutil)
        # Generate independent standard Laplace random variables
        Z <- matrix(rlaplace(N * P, m = 0, s = 1/sqrt(2)), nrow = N, ncol = P)
        
        Z <- Z / sqrt(2)
        
        L <- chol(sigma_0)
        data <- Z %*% L
    } else if(dist_type == "lognormal") { # Log-normal
        sdlog <- 1
        meanlog <- 0
        theoretical_mean <- exp(meanlog + sdlog^2/2)
        
        Z <- matrix(rlnorm(N * P, meanlog = meanlog, sdlog = sdlog), nrow = N, ncol = P)
        
        Z <- Z - theoretical_mean
        
        Z <- scale(Z, center = FALSE, scale = apply(Z, 2, sd))
        L <- chol(sigma_0)
        data <- Z %*% L
        }


    # sample covariance matrix
    S <- cov(data)
    diag(adj_mat_True) <- 1 
    total_edge_num <- sum(adj_mat_True) - P

    return(list(theta_0=theta_0, sigma_0=sigma_0, data=data, S=S, adj_mat_True=adj_mat_True, total_edge_num=total_edge_num, gtype=g))
}

args <- commandArgs(trailingOnly = TRUE)

# Parse the command line argument
nt <- args[1]
nt <- as.integer(nt)
cat("random seed set to", nt, "\n")

# get args for graph type and P
p_val <- Sys.getenv('P_VAL')
n_val <- Sys.getenv('N_VAL')
g_val <- Sys.getenv('G_VAL')
dist_type <- Sys.getenv('DIST_TYPE')
dist.Gtype <- strsplit(g_val,",")[[1]]
dist.p <- as.integer(strsplit(p_val, ",")[[1]])
n <- max(as.integer(strsplit(n_val, ",")[[1]]))

if(n <= 5000){
    n <- 5000
}

cat(dist.p,"\n",dist.Gtype, "\n")
cat("Generating", n, "sample data............\n")

# save dir
cur_dir <- getwd()
save_main_dir <- paste0(cur_dir, "/Data/")
# check if the directory exists
if(!dir.exists(save_main_dir)){
    dir.create(save_main_dir)
}

for(i in 1:length(dist.Gtype)){
    graph <- dist.Gtype[i]
    save_dir <- paste0(save_main_dir,graph,"_")
    for(p_idx in 1:length(dist.p)){
        P <- as.integer(dist.p[p_idx])
        G <- generator(N = n, P = P, g = graph, ss=nt, dist_type=dist_type)
        save_dir_p = paste0(save_dir, P, "_", nt, ".rds")
        cat("Saving dist for ", graph, " with dim = ", P, "in ", save_dir_p, "............")
        saveRDS(G, file = save_dir_p)
        cat("DONE\n")
    }
}


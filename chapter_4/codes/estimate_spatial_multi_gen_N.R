estimate_spatial_multi_gen_N <- function(Y, X, W, method = "s2sls",
                                         ind_beta = matrix(T, ncol(X), ncol(Y)),
                                         ind_RHO = matrix(T, ncol(Y), ncol(Y)),
                                         ind_GAMMA = matrix(F, ncol(Y), ncol(Y)), 
                                         compute_sd = F, instru = T) {
  
  # verification
  library("Matrix")
  stopifnot(method %in% c("s2sls", "s3sls"))
  
  # initialization
  M <- ncol(Y) # the dimension of Y
  n <- nrow(Y) # number of observations
  K <- ncol(X) # the size of X
  # check if the 1st column is the constant
  #if (K > 1) {
  #  if (! all(X[, 1] == 1)) {
  #    cat("Warning: Constant was not given in X")
 #     #X <- cbind(1, X)
 #     #K <- K + 1
#    }
#  }
  
  res_beta <- sd_beta <- matrix(0, K, M)
  hat_Uk <- matrix(0, n, M)
  RHO <- sd_RHO <-matrix(0, M, M)
  GAMMA <- sd_GAMMA <- matrix(0, M, M) 
  
  # verification
  stopifnot(nrow(X) == n)
  stopifnot(nrow(W) == n, ncol(W) == n)

  
  # initialisation
  if (instru) {
    W_X <- W %*% X[, -1]
    W2_X <- W %*% W_X
    H_n <- cbind(X, W_X, W2_X)
  } else {
    W_X <- W %*% X
    H_n <- cbind(X, W_X)    
  }
  W_Y <- W %*% Y
  
  P_H <- H_n %*% chol2inv(qr(H_n)$qr) %*% t(H_n)
  Z_tilde <- cbind(X, P_H %*% W_Y, P_H %*% Y)
  Z <- cbind(X, W_Y, Y)
  
    for (k in 1:M) {
      ind_beta_l <- ind_beta[, k]
      p_l <- length(which(ind_beta_l))
      ind_rho_l <- ind_RHO[k, ]
      rho_l <- length(which(ind_rho_l))
      ind_GAMMA_l <- ind_GAMMA[, k]
      GAMMA_l <- length(which(ind_GAMMA_l))
      ind_l <- c(ind_beta_l, ind_rho_l, ind_GAMMA_l) 
      
      Z_l <- Z[, ind_l]
      Z_tilde_l <- Z_tilde[, ind_l]
      
      if (k == 1) { # initialization
        Z_tilde_tot_l <- Z_tilde_l
        Z_tot_l <- Z_l
        } else { # we fill the vectorized matricial formula (3) in the article
          Z_tilde_tot_l <- Matrix::bdiag(Z_tilde_tot_l, Z_tilde_l)
          Z_tot_l <- Matrix::bdiag(Z_tot_l, Z_l)
        }
    
      cste <- chol2inv(qr(Z_tilde_l)$qr) %*% t(Z_tilde_l)
      res_k <- cste %*% Y[, k] 
      res_beta[ind_beta_l, k] <- res_k[1:p_l]
      RHO[k, ind_rho_l] <- res_k[(p_l + 1):(p_l + rho_l)]
      GAMMA[k, ind_GAMMA_l] <- res_k[(p_l + rho_l + 1):(p_l + rho_l + GAMMA_l)]
      hat_Uk[, k] <-  Y[, k] - Y %*% GAMMA[k, ] - X %*% res_beta[, k] - W_Y %*% RHO[k, ]
    }
  
  SIGMA <- t(hat_Uk) %*% hat_Uk / (n - K)
  
  if (method == "s3sls") {
    # if not, we start the s3sls estimation (formula (3))
    Sigtot <- kronecker(solve(SIGMA), diag(n))
    cste <- solve(t(Z_tilde_tot_l) %*% Sigtot %*% Z_tot_l) 
    res_k_3 <- cste %*% t(Z_tilde_tot_l) %*% Sigtot %*% as.vector(Y) 
    
    increment <- 0
    for (k in 1:M) {
      ind_beta_l <- ind_beta[, k]
      p_l <- length(which(ind_beta_l))
      ind_rho_l <- ind_RHO[k, ]
      rho_l <- length(which(ind_rho_l))
      ind_GAMMA_l <- ind_GAMMA[, k]
      GAMMA_l <- length(which(ind_GAMMA_l))
      ind_l <- c(ind_beta_l, ind_rho_l, ind_GAMMA_l) 
      
      res_beta[ind_beta_l, k] <- res_k_3[increment + (1:p_l)]
      RHO[k, ind_rho_l] <- res_k_3[increment + ((p_l + 1):(p_l + rho_l))]
      GAMMA[k, ind_GAMMA_l] <- res_k[increment + ((p_l + rho_l + 1):(p_l + rho_l + GAMMA_l))]
      
      increment <- increment + p_l + rho_l + GAMMA_l
      hat_Uk[, k] <-  Y[, k] - Y %*% GAMMA[k, ] - X %*% res_beta[, k] - W_Y %*% RHO[k, ]
    }
    
    SIGMA <- t(hat_Uk) %*% hat_Uk / (n - K)
  }
  
  if (compute_sd) {
    # compute the standard deviation of parameters
    Sigtot <- kronecker(solve(SIGMA), Diagonal(n))
    var_beta <- solve(crossprod(Z_tilde_tot_l, Sigtot) %*% Z_tot_l) # variance of beta hat
    std <- sqrt(diag(var_beta))
  
    increment <- 0
    for (k in 1:M) {
      ind_beta_l <- ind_beta[, k]
      p_l <- length(which(ind_beta_l))
      ind_rho_l <- ind_RHO[k, ]
      rho_l <- length(which(ind_rho_l))
      ind_GAMMA_l <- ind_GAMMA[, k]
      GAMMA_l <- length(which(ind_GAMMA_l))

      sd_beta[ind_beta_l, k] <- std[increment + (1:p_l)]
      sd_RHO[k, ind_rho_l] <- std[increment + (p_l + 1):(p_l + rho_l)]
      sd_GAMMA[k, ind_GAMMA_l] <- std[increment + (p_l + rho_l + 1):(p_l + rho_l + GAMMA_l)]
      increment <- increment + (p_l + rho_l + GAMMA_l)
    }
    
    return(list(est_beta = res_beta,
              est_RHO = RHO,
              est_GAMMA = GAMMA,
              est_SIGMA = SIGMA,
              sd_beta= sd_beta,
              sd_RHO = sd_RHO,
              sd_GAMMA = sd_GAMMA))
  } else {
    return(list(est_beta = res_beta,
                est_RHO = RHO,
                est_GAMMA = GAMMA,
                est_SIGMA = SIGMA
                ))
  }
}
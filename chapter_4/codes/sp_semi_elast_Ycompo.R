sp_semi_elast_Ycompo <- function(X, b_star, index_b, W, RHO, 
                                 method = "exact", type_x = "classic", V_x = NULL,
                                 delta = 10^(-8), print_summary = T, 
                                 V_y = ilrBase(ncol(index_b) + 1)) {
  
  # check 
  stopifnot(method %in% c("exact", "appro"))
  stopifnot(type_x %in% c("classic", "compo"))
  if (type_x == "compo")
    stopifnot(!is.null(V_x))
  
  # initialization
  n <- nrow(X)
  D <- ncol(b_star) + 1
  res <- array(0, dim = c(n, n, D))
  W_Mat <- as(W, "Matrix")
  A_w <- Diagonal(n * (D - 1)) - kronecker(RHO, W_Mat)
  A_w_inv <- solve(A_w)
  E_Y_ilr <- A_w_inv %*% as.vector(X %*% b_star)
  E_Y <- as(ilrInv(matrix(E_Y_ilr, ncol = D - 1), V = V_y), "matrix")
  
  if (method == "exact") {
    for (i in 1:n) {
      z <- E_Y[i, ]
      w <- (diag(D) - rep(1, D) %*% t(z)) %*% V_y # ilrBase(D = D)
      for (j in 1:n) {
        A_w_ij <- A_w_inv[c(i, i + n), c(j, j + n)]  
        temp <- w %*% A_w_ij %*% b_star[index_b, ]
        if (type_x == "compo")
          temp <- temp %*% V_x
        
        res[j, i, 1] <- temp[1]
        res[j, i, 2] <- temp[2]
        res[j, i, 3] <- temp[3]
      }
    } 
  } else {
    for (i in 1:n) {
      my_x_temp <- X
      my_x_temp[i, index_b] <- my_x_temp[i, index_b] + delta
      beta_x <- as.vector(my_x_temp %*% b_star)
      pred_y <- solve(A_w, beta_x)
      pred_y_simp <- as(ilrInv(matrix(pred_y, ncol = D - 1), V = V_y), "matrix")
      res[, i, 1] <- (pred_y_simp[, 1] - E_Y[, 1]) / (delta * E_Y[, 1])
      res[, i, 2] <- (pred_y_simp[, 2] - E_Y[, 2]) / (delta * E_Y[, 2])
      res[, i, 3] <- (pred_y_simp[, 3] - E_Y[, 3]) / (delta * E_Y[, 3])
    }
  }
  
  spatial_direct <- apply(res, 3, function(x) diag(x))
  spatial_total <- apply(res, c(1, 3), sum)
  spatial_indirect <- spatial_total - spatial_direct

  # if (print_summary) {
  #   cat("Direct impact: \n", 
  #     apply(E_Y * spatial_direct, 2, sum) / apply(E_Y, 2, sum), "\n",
  #     "Indirect impact: \n",
  #     apply(E_Y * spatial_indirect, 2, sum) / apply(E_Y, 2, sum), "\n",
  #     "Total impact: \n",
  #     apply(E_Y * spatial_total, 2, sum) / apply(E_Y, 2, sum), "\n"
  #     )
  # }
  if (print_summary) {
    cat("Direct impact: \n", 
        apply(spatial_direct, 2, sum) / n , "\n",
        "Indirect impact: \n",
        apply(spatial_indirect, 2, sum) / n, "\n",
        "Total impact: \n",
        apply(spatial_total, 2, sum) / n, "\n"
    )
  }                     
  return(res)
}

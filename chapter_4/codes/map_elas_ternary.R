map_elas_ternary <- function(df_coords, 
                             X, b_star, index_b, W, RHO, delta = 10^(-8),
                             type = 1, all_impacts = T, add = T, 
                             row_to_plot = 1:nrow(df_coords), scale_ternary = F,
                             ternary_unit = F, lwd = 1,
                             range, cex.pts = 0.5, size = 1, length = 0.1, cex.legend = 0.5,
                             legend.ternary = c("", "", ""),
                             names_ind = "", V_y = ilrBase(ncol(b_star) + 1)) {
  
  # initialisation 
  n <- nrow(df_coords)
  D <- ncol(b_star) + 1
  res <- array(0, dim = c(n, n, D))
  W_Mat <- as(W, "Matrix")
  A_w <- Diagonal(n * (D - 1)) - kronecker(RHO, W_Mat)
  A_w_inv <- solve(A_w)
  E_Y_ilr <- A_w_inv %*% as.vector(X %*% b_star)
  E_Y <- as(ilrInv(matrix(E_Y_ilr, ncol = D - 1), V = V_y), "matrix")
  
  # compute the semi-elasticities
  elast_sp <- sp_semi_elast_Ycompo(X, b_star, index_b = index_b, method = "appro",
                                   delta = delta, W = W, 
                                   RHO = RHO,
                                   print_summary = F, V_y = V_y)
  #elast_sp <- sp_semi_elast_Ycompo(X, b_star, index_b = index_b, method = "exact",
  #                                 W = W, 
  #                                 RHO = RHO,
  #                                 print_summary = F, V_y = V_y)
  # the prediction 
   pred_y_simp <- E_Y * (1 + delta * cbind(apply(elast_sp[ , , 1], 1, sum),
                                          apply(elast_sp[ , , 2], 1, sum),
                                          apply(elast_sp[ , , 3], 1, sum)))
   # equivalent to 
    X_temp <- X
    X_temp[, index_b] <- X_temp[, index_b] + delta
    E_Y_ilr_delta <- A_w_inv %*% as.vector(X_temp %*% b_star)
    pred_y_simp <- as(ilrInv(matrix(E_Y_ilr_delta, ncol = D - 1), V = V_y), 
                      "matrix")
  
  
  # the dimension of the triangle: s1 is the basis, s2 is the height
  s1 <- size * range / 10
  s2 <- s1 * sqrt(3) / 2
  
  if (!add) {
    if (ternary_unit) { 
      A <- c(0, 0)
      B <- c(1, 0)
      C <- c(0.5, sqrt(3) / 2)
      
      plot(c(A[1], B[1], C[1]), c(A[2], B[2], C[2]), 
           type = "n", xaxt = "n", yaxt = "n", ylim = c(-0.04, 0.9),
           xlab = "", ylab = "", frame = F, asp = 1)
      
      } else {
        plot(df_coords[row_to_plot, 1], df_coords[row_to_plot, 2], 
         type = "n", xaxt = "n", yaxt = "n", 
         xlab = "", ylab = "", frame = F, asp = 1)
      }
  }
  
  for (k in row_to_plot) {

    # define the individual triangle 
    if (!ternary_unit) { 
      x <- df_coords[k, 1]
      y <- df_coords[k, 2]
      A <- c(x - s1 / 2, y - s2 / 3)
      B <- c(x + s1 / 2, y - s2 / 3)
      C <- c(x, y + 2 * s2 / 3)

    } else {
      A <- c(0, 0)
      B <- c(1, 0)
      C <- c(0.5, sqrt(3) / 2)
    }

    lines(c(A[1], B[1]), c(A[2], B[2]), col = rgb(0.5, 0.5, 0.5)) # "#E16A86") # r
    lines(c(A[1], C[1]), c(A[2], C[2]), col = rgb(0.5, 0.5, 0.5)) # "#50A315") # g
    lines(c(B[1], C[1]), c(B[2], C[2]), col = rgb(0.5, 0.5, 0.5)) # "#009ADE") # b
    
    if (scale_ternary) {
      base_trian <- (B[1] - A[1]) 
      hauteur <- C[2] - A[2]
      # XR lines
      for (j in 0:5) {
        if(j != 0)
        lines(c(A[1] + base_trian * j / 10, B[1] - base_trian * j / 10), 
              c(A[2] + j / 5 * hauteur, A[2] + j / 5 * hauteur), 
              lty = 2, lwd = 0.7)
        text(B[1] - base_trian * j / 10, A[2] + j / 5 * hauteur, 2 * j / 10,
             pos = 4, cex = cex.legend)
      }
      
      # Left lines
      for (j in 0:5) {
        if(j != 0)
        lines(c(A[1] + j / 5 * base_trian, A[1] + base_trian / 2 + j * 1/10 * base_trian), 
              c(A[2], C[2] - j / 5 * hauteur), lty = 2, lwd = 0.7)
        text(A[1] + j / 5 * base_trian, A[2], 2 * j / 10, pos = 1, cex = cex.legend, srt = 70)
      }
      
      # Right lines 
      for (j in 0:5) {
        if(j != 5)
        lines(c(A[1] + base_trian * j * 1 / 10, A[1] + base_trian * j / 5), 
              c(A[2] + j / 5 * hauteur, A[2]), lty = 2, lwd = 0.7)
        text(A[1] + base_trian * j * 1 / 10, A[2] + j / 5 * hauteur, 1 - 2 * j / 10, 
             pos = 2, cex = cex.legend, srt = -50)
      }
    }
    
    # legend 
    hauteur <- C[2] - A[2]
    text(c(A[1]+0.02, B[1], C[1]), c(A[2] - hauteur / 40, B[2] - hauteur / 40, C[2] + hauteur / 40) + 0.005, 
         legend.ternary, pos = c(1, 1, 4), 
         cex = cex.legend, col = c("#C4A900", "#EE82E6", "#00C5C4"))
    
    # the adjusted values in the simplex
    E_Y_simplex_x <- E_Y[, 1] * A[1] + E_Y[, 2] * B[1] + E_Y[, 3] * C[1] 
    E_Y_simplex_y <- E_Y[, 1] * A[2] + E_Y[, 2] * B[2] + E_Y[, 3] * C[2]
    
    # the predicted values in the triangle
    Y_pred_x <- pred_y_simp[, 1] * A[1] + pred_y_simp[, 2] * B[1] + pred_y_simp[, 3] * C[1] 
    Y_pred_y <- pred_y_simp[, 1] * A[2] + pred_y_simp[, 2] * B[2] + pred_y_simp[, 3] * C[2]
    

    if (type == 1) {
      
      Y_simplex_x_start <- E_Y_simplex_x[k]
      Y_simplex_y_start <- E_Y_simplex_y[k]
      
      if (all_impacts) {
      
        for(l in 1:n) {
        
          pred_y_simp_ind <- c(
            elast_sp[k, l, 1] * (delta * E_Y[k, 1]) + E_Y[k, 1],
            elast_sp[k, l, 2] * (delta * E_Y[k, 2]) + E_Y[k, 2],
            elast_sp[k, l, 3] * (delta * E_Y[k, 3]) + E_Y[k, 3]
          )
        
          Y_simplex_x_end <- pred_y_simp_ind[1] * A[1] + pred_y_simp_ind[2] * B[1] +
            pred_y_simp_ind[3] * C[1]
          Y_simplex_y_end <- pred_y_simp_ind[1] * A[2] + pred_y_simp_ind[2] * B[2] +
            pred_y_simp_ind[3] * C[2]
        
          if(l != 1) {
            my_vec_x <- (Y_simplex_x_start - E_Y_simplex_x[k])
            my_vec_y <- (Y_simplex_y_start - E_Y_simplex_y[k])
            Y_simplex_x_end <- Y_simplex_x_end + my_vec_x
            Y_simplex_y_end <- Y_simplex_y_end + my_vec_y
          } else {
            text(Y_simplex_x_start, Y_simplex_y_start, names_ind, pos = 3, cex = 0.7)
          }
        
          arrows(Y_simplex_x_start, Y_simplex_y_start,
                 Y_simplex_x_end, Y_simplex_y_end,
                 col = ifelse(l == k, "#E16A86", "#009ADE"),
                 length = length, lwd = lwd)
        
          Y_simplex_x_start <- Y_simplex_x_end
          Y_simplex_y_start <- Y_simplex_y_end
        }
        #legend("topright", legend = c("direct", "indirect"), lty = 1,
        #       col = c("#E16A86", "#00AD9A"), cex = cex.legend)
      } else {
      
        # Direct impact
         pred_y_simp_ind <- c(
           elast_sp[k, k, 1] * (delta * E_Y[k, 1]) + E_Y[k, 1],
           elast_sp[k, k, 2] * (delta * E_Y[k, 2]) + E_Y[k, 2],
           elast_sp[k, k, 3] * (delta * E_Y[k, 3]) + E_Y[k, 3]
         )
         
         Y_simplex_x_end <- pred_y_simp_ind[1] * A[1] + pred_y_simp_ind[2] * B[1] +
           pred_y_simp_ind[3] * C[1]
         Y_simplex_y_end <- pred_y_simp_ind[1] * A[2] + pred_y_simp_ind[2] * B[2] +
           pred_y_simp_ind[3] * C[2]
         
         arrows(Y_simplex_x_start, Y_simplex_y_start,
                Y_simplex_x_end, Y_simplex_y_end,
                col = "#E16A86",
                length = length, lwd = lwd)
         
         # indirect impacts
         Y_simplex_x_start <- Y_simplex_x_end
         Y_simplex_y_start <- Y_simplex_y_end
         
         pred_y_simp_ind <- c(
           sum(elast_sp[k, -k, 1]) * (delta * E_Y[k, 1]) + E_Y[k, 1],
           sum(elast_sp[k, -k, 2]) * (delta * E_Y[k, 2]) + E_Y[k, 2],
           sum(elast_sp[k, -k, 3]) * (delta * E_Y[k, 3]) + E_Y[k, 3]
         )
         
         Y_simplex_x_end <- pred_y_simp_ind[1] * A[1] + pred_y_simp_ind[2] * B[1] +
           pred_y_simp_ind[3] * C[1]
         Y_simplex_y_end <- pred_y_simp_ind[1] * A[2] + pred_y_simp_ind[2] * B[2] +
           pred_y_simp_ind[3] * C[2]
         
         my_vec_x <- (Y_simplex_x_start - E_Y_simplex_x[k])
         my_vec_y <- (Y_simplex_y_start - E_Y_simplex_y[k])
         Y_simplex_x_end <- Y_simplex_x_end + my_vec_x
         Y_simplex_y_end <- Y_simplex_y_end + my_vec_y
         
  
         arrows(Y_simplex_x_start, Y_simplex_y_start,
                Y_simplex_x_end, Y_simplex_y_end,
                col = "#009ADE",
                length = length, lwd = lwd)
         
         #legend("topright", legend = c("total", "direct", "indirect"), lty = 1,
         #       col = c("#50A315", "#E16A86", "#009ADE"), cex = cex.legend)
         
         # total impact
         # points(Y_pred_x[k], Y_pred_y[k], pch = 16, cex = cex.pts, col = "#E16A86")
         arrows(E_Y_simplex_x[k], E_Y_simplex_y[k], Y_pred_x[k], Y_pred_y[k], 
                col = "#50A315", length = length, lwd = lwd)
      }
    } else {
      
      points(E_Y_simplex_x[-k], E_Y_simplex_y[-k], pch = 16, cex = cex.pts, 
             col = "#009ADE")
      for(l in c((1:n)[-k], k)) {
        
        pred_y_simp_ind <- c(
          elast_sp[l, k, 1] * (delta * E_Y[l, 1]) + E_Y[l, 1],
          elast_sp[l, k, 2] * (delta * E_Y[l, 2]) + E_Y[l, 2],
          elast_sp[l, k, 3] * (delta * E_Y[l, 3]) + E_Y[l, 3]
        )
        
        Y_simplex_x_end <- pred_y_simp_ind[1] * A[1] + pred_y_simp_ind[2] * B[1] +
          pred_y_simp_ind[3] * C[1]
        Y_simplex_y_end <- pred_y_simp_ind[1] * A[2] + pred_y_simp_ind[2] * B[2] +
          pred_y_simp_ind[3] * C[2]
        
        arrows(E_Y_simplex_x[l], E_Y_simplex_y[l], Y_simplex_x_end, Y_simplex_y_end, 
               col = ifelse(k == l, "#E16A86", "#009ADE"), length = length, lwd = lwd)
      }
      #arrows(E_Y_simplex_x[k], E_Y_simplex_y[k], Y_simplex_x_end, Y_simplex_y_end, 
      #       col = "#E16A86", length = length)
      #legend("topright", legend = c("direct", "indirect"), lty = 1,
      #       col = c("#E16A86", "#00AD9A"), cex = cex.legend)
    }
    

 
  }
  
}

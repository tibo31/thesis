map_elas_bar <- function(spatial_direct, spatial_indirect, spatial_total, coords_s,
                         index_to_plot = 1:nrow(spatial_direct),
                         q4 = colorspace::qualitative_hcl(ncol(spatial_direct), palette = "Dark 3"),
                         add = F, max_bar = 1, width_bar = 1, x.legend = "none", plot_scale = T, cex.legend = 1,
                         label_s = F, round.scale = 0, print_legend = T) {
  
  # verification
  # size of the vectors
 
  stopifnot(x.legend %in% c("none", "bottomright", "bottom", "bottomleft", 
            "left", "topleft", "top", "topright", "right", "center"))
  
  ############### Initialisation   
  # number of observation 
  n <- nrow(coords_s)
 
  # dimension of elasticities
  D <- ncol(spatial_direct)
  
  # check on coords_s
  stopifnot(ncol(coords_s) == 2)
  
  # check on elas
  stopifnot(nrow(coords_s) == nrow(spatial_direct))
  
  # coords must have rownames corresponding with S
  coords_s <- as.matrix(coords_s)
  
  # parameter to shift the destination barplot
  shift <- 1 / 50
    
  # site_s
  site_s <- rownames(coords_s)
  
  # define the coordinates of origin site: it corresponds to the coordinates
  # of s in S, slightly shifted
  shift_coords_x <- diff(range(coords_s[, 1])) * shift * width_bar
  shift_coords_x_2 <- shift_coords_x / 2
  
  # all eslasticities
  all_ela <- cbind(spatial_direct, spatial_indirect, spatial_total)
  
  x_direct <- matrix(coords_s[, 1], nrow = n, ncol = D) - shift_coords_x
  x_indirect <- matrix(coords_s[, 1], nrow = n, ncol = D) 
  x_total <- matrix(coords_s[, 1], nrow = n, ncol = D) + shift_coords_x
  
  # one bar per component 
  x_direct <- x_direct + matrix(seq(- shift_coords_x_2 / 2, shift_coords_x_2 / 2, length.out = D),
                                nrow = n, ncol = D, byrow = T)
  x_indirect <- x_indirect + matrix(seq(- shift_coords_x_2 / 2, shift_coords_x_2 / 2, length.out = D),
                                nrow = n, ncol = D, byrow = T)
  x_total <- x_total + matrix(seq(- shift_coords_x_2 / 2, shift_coords_x_2 / 2, length.out = D),
                                nrow = n, ncol = D, byrow = T)
  
  # width of the bar 
   shift_bar <- diff(seq(- shift_coords_x_2 / 2, shift_coords_x_2 / 2, length.out = D)[1:2]) / 2
    
  # maximum height for the bars 
  max_bar <- diff(range(coords_s[, 2])) * shift * max_bar
  bar_direct <-  max_bar * spatial_direct / max(abs(all_ela), na.rm = T) 
  bar_indirect <- max_bar * spatial_indirect / max(abs(all_ela), na.rm = T) 
  bar_total <- max_bar * spatial_total / max(abs(all_ela), na.rm = T) 
  
  if (!add) {  
   plot(coords_s[, 1], coords_s[, 2], type = "n", xaxt = "n", yaxt = "n", 
        xlab = "", ylab = "", frame = F, asp = 1)
  }

  # outflows bar
  for(k in index_to_plot) {
    # plot the horizontal lines
    values_hori <- round(c(max(abs(all_ela)), max(abs(all_ela)) / 2, 0, 
                           - max(abs(all_ela)) / 2, - max(abs(all_ela))), 
                         round.scale)
    lines(c(x_direct[k, 1] - shift_bar, x_total[k, D] + shift_bar), 
          coords_s[k, 2] + max_bar * c(values_hori[1], values_hori[1]) / max(abs(all_ela), na.rm = T), 
          lty = 1, col = "white")
    lines(c(x_direct[k, 1] - shift_bar, x_total[k, D] + shift_bar), 
          coords_s[k, 2] + max_bar * c(values_hori[2], values_hori[2]) / max(abs(all_ela), na.rm = T), 
          lty = 1, col = "white")
    lines(c(x_direct[k, 1] - shift_bar, x_total[k, D] + shift_bar), 
          coords_s[k, 2] + max_bar * c(values_hori[3], values_hori[3]) / max(abs(all_ela), na.rm = T), 
          lty = 1, col = "white")
    lines(c(x_direct[k, 1] - shift_bar, x_total[k, D] + shift_bar), 
          coords_s[k, 2] + max_bar * c(values_hori[4], values_hori[4]) / max(abs(all_ela), na.rm = T), 
          lty = 1, col = "white")
    lines(c(x_direct[k, 1] - shift_bar, x_total[k, D] + shift_bar), 
          coords_s[k, 2] + max_bar * c(values_hori[5], values_hori[5]) / max(abs(all_ela), na.rm = T), 
          lty = 1, col = "white")
   
    if (plot_scale) {
      #points(rep(x_direct[k, 1] - shift_bar, 5), 
      #       coords_s[k, 2] + round(c(max(abs(all_ela)), max(abs(all_ela)) / 2, 0, 
      #                                - max(abs(all_ela)) / 2, - max(abs(all_ela))), 
      #                              round.scale) * max_bar / max(abs(all_ela), na.rm = T),
      #       pch = "-", cex = 0.2)
      text(rep(x_direct[k, 1] - shift_bar, 5), 
           coords_s[k, 2] + round(c(max(abs(all_ela)), max(abs(all_ela)) / 2, 0, 
                                    - max(abs(all_ela)) / 2, - max(abs(all_ela))), 
                                  round.scale) * max_bar / max(abs(all_ela), na.rm = T),
           values_hori,cex = cex.legend, pos = 2)
    }
    for (j in 1:D) {
      # direct 
      bar_x <- cbind(x_direct[k, j] - shift_bar, x_direct[k, j] + shift_bar, 
                     x_direct[k, j] + shift_bar, x_direct[k, j] - shift_bar,
                     x_direct[k, j] - shift_bar)
      bar_y <- cbind(coords_s[k, 2], coords_s[k, 2], 
                      coords_s[k, 2] + bar_direct[k, j], 
                      coords_s[k, 2] + bar_direct[k, j],
                      coords_s[k, 2])
      polygon(bar_x, bar_y, col = q4[j])
      # indirect
      bar_x <- cbind(x_indirect[k, j] - shift_bar, x_indirect[k, j] + shift_bar, 
                     x_indirect[k, j] + shift_bar, x_indirect[k, j] - shift_bar,
                     x_indirect[k, j] - shift_bar)
      bar_y <- cbind(coords_s[k, 2], coords_s[k, 2], 
                     coords_s[k, 2] + bar_indirect[k, j], 
                     coords_s[k, 2] + bar_indirect[k, j],
                     coords_s[k, 2])
      polygon(bar_x, bar_y, col = q4[j])
      # total 
      bar_x <- cbind(x_total[k, j] - shift_bar, x_total[k, j] + shift_bar, 
                     x_total[k, j] + shift_bar, x_total[k, j] - shift_bar,
                     x_total[k, j] - shift_bar)
      bar_y <- cbind(coords_s[k, 2], coords_s[k, 2], 
                     coords_s[k, 2] + bar_total[k, j], 
                     coords_s[k, 2] + bar_total[k, j],
                     coords_s[k, 2])
      polygon(bar_x, bar_y, col = q4[j])
    }
   

  # print labels 
  if (label_s) {
   #bar_out_S <-  max_bar * outflows / max(c(inflows, outflows), na.rm = T) 
   #bar_in_S <- max_bar * inflows / max(c(inflows, outflows), na.rm = T) 
   text(coords_s[k, 1] , coords_s[k, 2] + max_bar, site_s[k], pos = 3, cex = cex.legend)
  }
  
  # plot the legend
  if (x.legend != "none" & print_legend) {
    text(coords_s[k, 1] - shift_coords_x, 
         coords_s[k, 2] - max_bar/1.25, "DE", srt = 0, cex = cex.legend)
    text(coords_s[k, 1], 
         coords_s[k, 2]  - max_bar/1.25, "IE", srt = 0, cex = cex.legend)
    text(coords_s[k, 1] + shift_coords_x, 
         coords_s[k, 2] - max_bar/1.25, "TE", srt = 0, cex = cex.legend)
    if (is.null(colnames(spatial_direct))) {
      names_D <- paste0("y", 1:D)
    } else {
      names_D <- colnames(spatial_direct) 
    }
  }
  }
  if(print_legend)
    legend(x.legend, legend = names_D, fill = q4,  box.lwd = -4, 
           cex = cex.legend)
}


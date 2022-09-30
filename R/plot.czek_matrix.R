## This file is part of RMaCzek

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


#'@title Produce a Czekanowski's Diagram
#'@description This is a function that can produce a Czekanowski's Diagram and present clustering findings.
#'@param x  a matrix with class czek_matrix.
#'@param type specifies if the graph should use color or symbols. The standard setting is symbols.
#'@param values specifies the color or the size of the symbols in the graph. The standard setting is a grey scale for a color graph and a vector with the values 2,1,0.5,0.25 and 0 for a graph with symbols.
#'@param plot_title specifies the main title in the graph.
#'@param tl.cex  Numeric, for the size of text label.
#'@param tl.offset Numeric, for text label.
#'@param tl.srt Numeric, for text label, string rotation in degrees.
#'@param pal The colour vector representing the clusters.
#'@param alpha  Factor modifying the opacity, alpha, typically in [0,1].
#'@param ps_power A power value to adjust point size.
#'@param col_size When type="col", the size of each point (maximum is 1).
#'@param cex.main Specify the size of the title text.
#'@param ... specifies further parameters that can be passed on to the plot function.
#'@export
#'@examples
#'# Set data ####
#'# Not Cluster
#'czek = czek_matrix(mtcars)
#'# Exact Clustering
#'czek_exact = czek_matrix(mtcars, order = "GW", cluster = TRUE, num_cluster = 2, min.size = 2)
#'# Fuzzy Clustering
#'czek_fuzzy = czek_matrix(mtcars, order = "OLO", cluster = TRUE, num_cluster = 2,
#'cluster_type = "fuzzy", min.size = 2, scale_bandwidth = 0.2)
#'
#'# Standard plot ############
#'plot(czek_exact)
#'plot.czek_matrix(czek_fuzzy)
#'
#'# Edit diagram title
#'plot(czek, plot_title = "mtcars", cex.main = 2)
#'
#'# Change point size ############
#'# Specify values
#'plot(czek, values = c(1, 0.8, 0.5, 0.2, 0))
#'plot(czek, values = grDevices::colorRampPalette(c("black", "red", "white"))(5))
#'
#'# set point size for 'symbols' type by setting power value
#'plot(czek, type = "symbols", ps_power = 1)
#'
#'# set point size for 'col' type
#'plot(czek, type = "col", col_size = 0.6)
#'
#'# Specify type ############
#'plot(czek, type = "symbols")
#'plot(czek, type = "col")
#'
#'# Edit cluster ############
#'# Edit colors
#'plot(czek_exact, pal = c("red", "blue"))
#'# Edit opacity
#'plot(czek_exact, alpha = 0.5)
#'
plot.czek_matrix = function(x, values = NULL, type = "symbols", plot_title = "Czekanowski's diagram",
                            tl.cex = 1, tl.offset = 0.4, tl.srt = 90,
                            pal = brewer.pal(n = 8, name = "Dark2"), alpha = 0.3, ps_power = 0.6,
                            col_size = 1, cex.main = 1, ...){

  oldpar = par(mar = c(0, 0, 4, 0))
  on.exit(par(oldpar))

  n_classes <- attr(x, "n_classes")
  levels <- attr(x, "levels")
  partition_boundaries <- attr(x, "partition_boundaries")
  new_order <- attr(x, "order")
  cluster <- attr(x, "cluster")
  cluster_boundary <- attr(x, "cluster_boundary")
  names = colnames(x)[new_order]
  x <- x[new_order, new_order]

  if (class(values) %in% c("numeric", "character") & length(values) == n_classes) {
    values <- values
  }else if (type == "symbols") {
    values <- rep(0, n_classes)
    values[1] <- 1
    for (i in 2:(n_classes - 1)) {
      values[i] <- values[i - 1]/2
    }
    values[n_classes] <- 0
  }else if (type == "col") {
    values <- round(seq(0, 100, length.out = n_classes))
    values <- paste("gray", values, sep = "")
  }else stop("type should be either 'col' or 'symbols'")

  n <- nrow(x)
  plot_values <- values[x]
  plot_y <- rep(n:1, n)
  plot_x <- rep(1:n, rep(n, n))

  graphics::plot.new()
  xlabwidth = max(graphics::strwidth(names, cex = tl.cex))
  ylabwidth = max(graphics::strwidth(names, cex = tl.cex))
  laboffset = graphics::strwidth("W", cex = tl.cex)

  for (i in 1:50) {
    xlim = c(1 - 0.5 - laboffset - xlabwidth,
             n + 0.5 + xlabwidth * abs(cos(tl.srt * pi/180)))
    ylim = c(0 - 0.5 - laboffset - ylabwidth * abs(sin(tl.srt * pi/180)) - tl.offset,
             n + 0.5 + laboffset)

    graphics::plot.window(xlim, ylim, asp = 1, xaxs = "i", yaxs = "i")

    x.tmp = max(graphics::strwidth(names, cex = tl.cex))
    y.tmp = max(graphics::strwidth(names, cex = tl.cex))
    laboffset.tmp = graphics::strwidth("W", cex = tl.cex) * tl.offset
    if (max(x.tmp - xlabwidth, y.tmp - ylabwidth, laboffset.tmp - laboffset) < 0.001) {
      break
    }
    xlabwidth = x.tmp
    ylabwidth = y.tmp
    laboffset = laboffset.tmp
    if (i == 50) {
      warning(c("Not been able to calculate text margin, ",
                "please try again with a clean new empty window using ",
                "{plot.new(); dev.off()} or reduce tl.cex"))
    }
  }

  if (.Platform$OS.type == "windows") {
    grDevices::windows.options(width = 7, height = 7 * diff(ylim)/diff(xlim))
  }

  plot_y <- rep(n:1, n)
  plot_x <- rep(1:n, rep(n, n))

  if (is.character(plot_values)) {
    graphics::symbols(plot_x, plot_y, add = TRUE, inches = FALSE, squares = rep(col_size, length(plot_x)), fg = plot_values, bg = plot_values)
  }else if (is.numeric(plot_values)) {
    graphics::symbols(plot_x, plot_y, add = TRUE, inches = FALSE, circles = plot_values^ps_power/2, bg = "black")

    # hide the zero points
    ind.p = which(plot_values == 0)
    if(length(ind.p) != 0){
      graphics::symbols(plot_x[ind.p], plot_y[ind.p],
              inches = FALSE, squares = rep(1, length(plot_x[ind.p])), fg = "white", bg = "white", add = TRUE)
    }
  }

  graphics::rect(0.5, 0.5, n + 0.5, n + 0.5, border = "black")

  cex = seq(0.7, 1, by = 0.1)
  height = sapply(cex, function(cex) graphics::strheight("A", cex = cex)) * n
  short_aes = min(diff(xlim), diff(ylim))
  tl.cex = cex[which.min(abs(height - short_aes))]

  # x and y labels
  graphics::axis(1, at = 1:n, tick = FALSE, labels = names, las = 2, cex.axis = tl.cex, pos = tl.offset)
  graphics::axis(2, at = n:1, tick = FALSE, labels = names, las = 1, cex.axis = tl.cex, pos = tl.offset)

  graphics::title(main = plot_title, cex.main = cex.main)

  if(!is.null(cluster)){
    num_cluster = nrow(cluster_boundary)
    for (i in 1:num_cluster) {
      # rect(cluster_boundary[i,1], cluster_boundary[i,2], cluster_boundary[i,3],cluster_boundary[i,4], col= adjustcolor(pal[i], alpha.f = alpha), border = FALSE)
      graphics::rect(cluster_boundary[i,1], cluster_boundary[i,2], cluster_boundary[i,3],cluster_boundary[i,4], col= grDevices::adjustcolor(pal[i], alpha.f = alpha), border = FALSE)
    }
  }

  # for(i in 1:n){
  #   # symbols(i, n-i+1, add = TRUE, inches = FALSE, circles = 0.5, bg = col[4-true_label[i]], fg = col[4-true_label[i]])
  #   symbols(i, n-i+1, inches = FALSE, squares = 1, fg = col[4-true_label[i]], bg = col[4-true_label[i]], add = TRUE)
  # }
}


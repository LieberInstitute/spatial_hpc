var.n <- function(array) {
  variance <- var(array)*(length(array)-1)/length(array)
  return(variance)
}


sd.n <- function(array) {
  return(sqrt(var.n(array)))
}

calculate.adj.matrix <- function(x, y, x_pixel = NULL, y_pixel = NULL, image = NULL,
                                 beta = 49, alpha = 1, histology = TRUE){
  # x, y, x_pixel, y_pixel are lists
  if (histology == TRUE){
    stopifnot(!is.null(x_pixel) & !is.null(y_pixel) & !is.null(image))
    stopifnot((length(x) == length(x_pixel)) & (length(y) == length(y_pixel)))
    cat("Calculating adj matrix using histology image...")
    # beta to control the range of neighborhood when calculate grey value for one spot
    # alpha to control the color scale
    beta_half <- round(beta/2)
    g <- NULL
    for (i in 1:length(x_pixel)){
      max_x <- dim(image)[1]
      max_y <- dim(image)[2]
      # need to round x_pixel and y_pixel for calculating adj matrix for low or high resolution image
      # pixel starts with 0
      first1 = as.character(max(1,x_pixel[i]-beta_half+1))
      first2 = as.character(min(max_x, x_pixel[i]+beta_half+1))
      second1 = as.character(max(1,y_pixel[i]-beta_half+1))
      second2 = as.character(min(max_y, y_pixel[i]+beta_half+1))
      a = as.numeric(first1)
      b = as.numeric(first2)
      c = as.numeric(second1)
      d = as.numeric(second2)
      nbs <- image[a:b, c:d,]
      # create average rgb for each block
      g <- rbind(g, apply(nbs,2,mean))
    }
    
    var.r <- var.n(g[,1])
    var.g <- var.n(g[,2])
    var.b <- var.n(g[,3])
    print(paste("Var of r, g, b = ", var.r, var.g, var.b))
    # raw score
    z0 <- (g[,1]*var.r+g[,2]*var.g+g[,3]*var.b)/(var.r + var.g + var.b)
    # scaled z score
    z1 <- (z0 - mean(z0))/sd.n(z0)
    z_scale <- max(sd.n(x), sd.n(y))*alpha
    z <- z1*z_scale
    print(paste("Var of x, y, z = ", var.n(x), var.n(y), var.n(z)))
    mat <- as.data.frame(cbind(x, y, z))
  }else {
    cat("Calculating adj matrix using xy only...")
    mat <- as.data.frame(cbind(x, y))
  }
  adj <- as.matrix(dist(mat, method = "euclidean", diag = TRUE, upper = TRUE))
  return(adj)
}

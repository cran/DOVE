# Internal function: plot the curves of VE_a and VE_h 
# knots and tau are provided in days
#' @importFrom graphics par plot lines axis legend
#' @importFrom grDevices dev.new
#' @importFrom stats stepfun
#' @importFrom utils tail
# 6/30/2021 -- extracted from plot.doveObj of v1.6 to
# allow for alternative implementation for dove2()
# 7/30/2021 -- constrain VE_h to be 0 at day 0 
# and non-decreasing at the right tail (lines 92-116)
.VEplot_type1 <- function(x, bandwidth = NULL, ...) {

  ve <- x$vaccine[[ 1L ]]
  
  if (is.null(x = bandwidth)) bandwidth <- 0.3
  
  if (!is.numeric(x = bandwidth)) {
    stop("bandwidth must be numeric", call. = FALSE)
  }
  
  # plot vaccine efficacy in percent as a step function
  
  opts20 <- seq(-200, 200, 20)
  opts10 <- seq(-200, 200, 10)

  ymin <- min(floor(x = floor(x = min(ve[,c(2L,4L,5L)]*100.0)) * 1.1),0.0)
  tst <- which(x = opts20 > ymin)[1L]
  if (tst > 1L) {
    ymin20 <- opts20[which(x = opts20 > ymin)[1L] - 1L]
    ymin10 <- opts10[which(x = opts10 > ymin)[1L] - 1L]
  } else {
    ymin20 <- ymin
    ymin10 <- ymin
  }

  ymax <- ceiling(x = ceiling(x = max(ve[,c(2L,4L,5L)]*100.0)) * 1.1)
  ymax20 <- min(opts20[which(x = opts20 > ymax)[1L]], 100.0)
  ymax10 <- min(opts10[which(x = opts10 > ymax)[1L]], 100.0)

  if (ymin20 == ymin10 && ymax20 == ymax10) {
    ymin <- ymin20
    ymax <- ymax20
    dy <- 20L
  } else {
    ymin <- ymin10
    ymax <- ymax10
    dy <- 10L
  }
  
  xmax <- ceiling(x = max(ve[,1L]))
  
  stepMain <- stats::stepfun(x =  ve[-1L,1L], y = ve[,2L]*100.0)
  
  stats::plot.stepfun(x = stepMain,
                      xlab = "Time Since Vaccination (in Days)", 
                      ylab = "Vaccine Efficiency in Reducing Attack Rate (%)",
                      xlim = c(0.0, xmax), ylim = c(ymin, ymax),
                      yaxt = "n", do.points = FALSE, main = "")
  
  graphics::axis(side = 2L, 
                 at = seq(from = ymin, to = ymax, by = dy), 
                 labels = paste0(seq(from = ymin, to = ymax, by = dy), "%"), 
                 cex = 0.8)
  
  stepLow <- stats::stepfun(x = ve[-1L,1L], y = ve[,4L]*100.0)
  
  graphics::lines(x = stepLow, col = 3L, do.points = FALSE)
  
  stepHigh <- stats::stepfun(x = ve[-1L,1L], y = ve[,5L]*100.0)
  
  graphics::lines(x = stepHigh, col = 3L, do.points = FALSE)
  
  graphics::legend(x = "bottom", 
                   legend = c("VE_a", "95% CI"), 
                   lty = c(1,1),
                   col = c(1L,3L), 
                   bg = "gray95")
  
  grDevices::dev.new()
  
  # reduction in hazard ratio as a percent
  
  testPoints <- c(0.0, attr(x,"tau")*{1L:100L}/100.0)
  
  nbw <- length(x = bandwidth)
  
  aa <- matrix(data = 0.0, nrow = length(x = testPoints), ncol = nbw)
  
  for (i in 1L:nbw) {
    
    bw <- bandwidth[i]*
      {max(ve[,1L]) - min(ve[,1L])} / 
      {attr(x,"dSum")^(1.0/5.0)}

    dmatrix <- outer(X = testPoints, 
                     Y = ve[,1L], 
                     FUN = "-") / bw

    kmatrix <- exp(x = -dmatrix^2L / 2.0) / 
      bw / sqrt(x = 2.0*pi)

    for (j in 1L:length(x = testPoints)) {

      weight <- sqrt(x = kmatrix[j,])

      # {npt x 2}
      Xj <- cbind(weight, weight*dmatrix[j,])
      # {np2}
      Yj <- weight*ve[,6L]
      
      tempj <- solve(a = crossprod(x = Xj), b = crossprod(x = Xj, y = Yj))

      aa[j,i] <- tempj[1L]*sum(kmatrix[j,])
    }

    if (any(aa[,i] <= 1e-12)) {
      ind <- which(x = testPoints < {28.0+1e-8})
      tst <- aa[ind,i] <= 1e-12
      aa[ind,i][tst] <- 1e-8
    }

  }

  # 7/30/21: Remove the estimates that are less than 28 days
  # and linearly extrapolate the estimates between day 0 and day 28
  # on the log scale
  ind <- which(x = testPoints < {28.0+1e-8})
  if (length(x = ind) > 1L) {
    tmp <-  utils::tail(x = ind, n = 1L)

    for (i in 1L:nbw) {
      slope <- log(x = aa[tmp,i]) / testPoints[tmp]
      aa[ind,i] <- exp(x = slope*testPoints[ind])
    }
  } else if (length(x = ind) == 1L) {
    aa[ind,] <- 1.0
  }

  # 7/30/21: make the right tail non-increasing by extrapolating 
  # the smallest value to the end.
  for (i in 1L:nbw) {
    tempdiff <- diff(x = aa[,i])
    tempind <- which(x = tempdiff >= 1e-8)
    if (length(x = tempind) > 0L) {
      right.tail.ind <- utils::tail(x = tempind, n = 1L)
      if (right.tail.ind < length(x = testPoints)) {
        aa[right.tail.ind:length(x = testPoints),i] <- aa[right.tail.ind,i]
      }
    }
  }
  ymin <- min(floor(x = floor(x = min({1.0-aa}*100.0)) * 1.1),0.0)
  ymin20 <- opts20[which(x = opts20 > ymin)[1L] - 1L]
  ymin10 <- opts10[which(x = opts10 > ymin)[1L] - 1L]

  ymax <- ceiling(x = ceiling(x = max({1.0-aa}*100.0)) * 1.1)
  ymax20 <- min(opts20[which(x = opts20 > ymax)[1L]], 100.0)
  ymax10 <- min(opts10[which(x = opts10 > ymax)[1L]], 100.0)

  if (ymin20 == ymin10 && ymax20 == ymax10) {
    ymin <- ymin20
    ymax <- ymax20
    dy <- 20L
  } else {
    ymin <- ymin10
    ymax <- ymax10
    dy <- 10L
  }
  
  
  graphics::plot(x = testPoints, y = {1.0-aa[,1L]}*100,
                 xlab = "Time Since Vaccination (in Days)", 
                 ylab = "Vaccine Efficiency in Reducing Hazard Rate (%)", 
                 type = 'l',
                 yaxt="n", main="", ylim = c(ymin, ymax))
  
  i <- 2L
  while (i <= nbw) {
    graphics::lines(x = testPoints, y = {1.0-aa[,i]}*100, col = i)
    i <- i + 1L
  }
  
  graphics::axis(side = 2L, 
                 at = seq(from = ymin, to = ymax, by = dy), 
                 labels = paste0(seq(from = ymin, to = ymax, by = dy), "%"), 
                 cex = 0.8)
  
  graphics::legend(x = "bottom", 
                   legend = round(x = bandwidth, digits = 4), 
                   lty = rep(x = 1L, times = nbw),
                   col = 1L:nbw, 
                   bg = "gray95",
                   title = "bandwidth")
  
  return()
  
}

# Plot Analysis Results
#
# Plot the key results of the analysis.
#
# @param x A doveObj object. The value returned by dove().
#
# @param y Ignored.
#
# @param ... Ignored. 
#
# @param bandwidth A numeric vector object. A tuning parameter for the 
#    bandwidth used for kernel estimation of the vaccine efficacy in 
#    reducing the hazard rate; this input is ignored if plots=FALSE.
#
#
#
# @name plot
# @examples
# data(doveData)
#
# set.seed(1234)
#
# ind <- sample(1:nrow(x = doveData), 2500, FALSE)
#
# # NOTE: This sample size is chosen for example only -- larger data sets
# # should be used.
# # See the vignette for a full analysis of the doveData dataset
#
# res <- dove(formula = Surv(event.time, event.status) ~ priority + sex + 
#                        vaccine(entry.time, vaccine.status, vaccine.time),
#             data = doveData[ind,])
#
# plot(x = res, bandwidth = c(0.5,1.0))
# @method plot doveObj
# @export 
plot.doveObj <- function(x, y, ..., bandwidth = NULL) {

  if (is.null(x = bandwidth)) bandwidth <- 0.5
  if (!is.numeric(x = bandwidth)) {
    stop("bandwidth must be numeric", call. = FALSE)
  }

  # plot vaccine efficacy in percent as a step function

  ymin <- min(x$vaccine[[ 1L ]][,c(2L,4L,5L)]*100.0)
  ymin <- ymin - ymin %% 10L - 10.0
    
  xmax <- ceiling(x = max(x$vaccine[[ 1L ]][,1L]))

  stepMain <- stats::stepfun(x =  x$vaccine[[ 1L ]][-1L,1L], 
                             y = x$vaccine[[ 1L ]][,2L]*100.0)
    
  stats::plot.stepfun(x = stepMain,
                      xlab = "Time Since Vaccination (in Days)", 
                      ylab = "Vaccine Efficiency in Reducing Attack Rate (%)",
                      xlim = c(0.0, xmax), ylim = c(ymin, 100),
                      yaxt = "n", do.points = FALSE, main = "")
    
  graphics::axis(side = 2L, 
                 at = seq(from = ymin, to = 100, by = 10), 
                 labels = paste0(seq(from = ymin, to = 100, by = 10), "%"), 
                 cex = 0.8)
    
  stepLow <- stats::stepfun(x = x$vaccine[[ 1L ]][-1L,1L], 
                            y = x$vaccine[[ 1L ]][,4L]*100.0)

  graphics::lines(x = stepLow, col = 3L, do.points = FALSE)

  stepHigh <- stats::stepfun(x = x$vaccine[[ 1L ]][-1L,1L], 
                             y = x$vaccine[[ 1L ]][,5L]*100.0)

  graphics::lines(x = stepHigh, col = 3L, do.points = FALSE)
    
  graphics::legend(x = "bottomleft", 
                   legend = c("VE_a", "95% CI"), 
                   lty = c(1,1),
                   col = c(1L,3L), 
                   bg = "gray95")

  grDevices::dev.new()
    
  # reduction in hazard ratio as a percent

  testPoints <- c(0.0, attr(x,"tau")*{1L:100L}/100.0)

  nbw <- length(x = bandwidth)

  aa <- matrix(data = 0.0, 
               nrow = length(x = testPoints),  
               ncol = nbw)

  for (i in 1L:nbw) {

    bw <- bandwidth[i]*
          {max(x$vaccine[[ 1L ]][,1L]) - min(x$vaccine[[ 1L ]][,1L])} / 
          {attr(x,"dSum")^(1.0/5.0)}

    dmatrix <- outer(X = testPoints, 
                     Y = x$vaccine[[ 1L ]][,1L], 
                     FUN = "-") / bw

    kmatrix <- exp(x = -dmatrix^2L / 2.0) / 
               bw / sqrt(x = 2.0*pi)

    for (j in 1L:length(x = testPoints)) {
      # {npt x 2}
      Xj <- cbind(kmatrix[j,], kmatrix[j,]*dmatrix[j,])
      # {np2}
      Yj <- kmatrix[j,]*x$vaccine[[ 1L ]][,6L]

      tempj <- solve(a = crossprod(x = Xj), b = crossprod(x = Xj, y = Yj))

      aa[j,i] <- tempj[1L]*sum(kmatrix[j,])
    }
  }

  ymin <- min({1.0-aa}*100)
  ymin <- ymin - ymin %% 10L - 10.0

  graphics::plot(x = testPoints, y = {1.0-aa[,1L]}*100,
                 xlab = "Time Since Vaccination (in Days)", 
                 ylab = "Vaccine Efficiency in Reducing Hazard Rate (%)", 
                 type = 'l',
                 yaxt="n", main="", ylim = c(ymin, 100))

  i <- 2L
  while (i <= nbw) {
    graphics::lines(x = testPoints, y = {1.0-aa[,i]}*100, col = i)
    i <- i + 1L
  }
    
  graphics::axis(side = 2L, 
                 at = seq(from = ymin, to = 100, by = 10), 
                 labels = paste0(seq(from = ymin, to = 100, by = 10), "%"), 
                 cex = 0.8)

  if (nbw > 1L) {
    graphics::legend(x = "bottomleft", 
                     legend = round(x = bandwidth, digits = 4), 
                     lty = rep(x = 1L, times = nbw),
                     col = 1L:nbw, 
                     bg = "gray95")
  }

  return()


}


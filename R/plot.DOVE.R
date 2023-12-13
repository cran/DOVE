#' Plot Estimated Vaccine Efficacy
#'
#' Generates plots of the estimated vaccine efficacy in reducing attack rate, 
#'   the estimated vaccine efficacy in reducing the hazard rate, 
#'   and their 95\% confidence intervals.
#'
#' @param x An DOVE object. The value object returned by dove() and dove2().
#'
#' @param y ignored.
#'
#' @param ... ignored
#' 
#' @param bandwidth A numeric vector object. A tuning parameter for the 
#'    bandwidth used for kernel estimation of the vaccine efficacy in 
#'    reducing the hazard rate for objects returned by dove(). This
#'    input is ignored for objects returned by dove2().
#' 
#' @name plot
#' @method plot DOVE
#'
#' @returns No return value, called to produce graphical elements.
#'
#' @examples
#' \dontshow{
#'   RcppArmadillo::armadillo_throttle_cores(2)
#' }
#' data(doveData)
#'
#' set.seed(1234)
#'
#' smp <- sample(1:nrow(x = doveData), 500, FALSE)
#'
#' # NOTE: This sample size is chosen for example only -- larger data sets
#' # should be used.
#' # See the vignette for a full dove() analysis of the doveData dataset
#'
#' result <- dove(formula = Surv(event.time, event.status) ~ priority + sex + 
#'                          vaccine(entry.time, vaccine.status, vaccine.time),
#'                data = doveData[smp,])
#'
#' plot(x = result, bandwidth = c(0.5,1.0))
#'
#' set.seed(1234)
#' smp <- sample(1L:nrow(x = doveData), size = 2500L)
#' 
#' # NOTE: This sample size is chosen for example only -- larger data sets
#' # should be used.
#' # See the vignette for a full dove2() analysis of the doveData dataset
#'
#' # Fit the model with default settings
#' result <- dove2(formula = Surv(event.time, event.status) ~ priority + sex + 
#'                           vaccine(entry.time, vaccine.status, vaccine.time), 
#'                 data = doveData[smp,])
#' 
#' plot(x = result)
#'
#' @export 
#' @include VEplot_type1.R VEplot_type2.R
# 6/30/2021 modified from original to include new dove2() plotting
# capabilities
plot.DOVE <- function(x, y, ..., bandwidth = NULL) {

  if (attr(x = x, which = "type") == 3L) {

    message("plot() is not available for the provided analysis ",
            "(constantVE = TRUE)")

  } else if (attr(x = x, which = "type") == 2L) {
    
    # 8/12/21: change the argument of .VEplot_type2

    # .VEplot_type2(knots = attr(x = x, which = "knots"),
    #               tau = attr(x = x, which = "tau"),
    #               gamma = attr(x = x, which = "gamma"),
    #               covgamma = attr(x = x, which = "covgamma"))
    
    .VEplot_type2(VE_a = x$vaccine$VE_a, VE_h = x$vaccine$VE_h)

  } else if (attr(x = x, which = "type") == 1L) {

    .VEplot_type1(x = x, bandwidth = bandwidth)

  }

}


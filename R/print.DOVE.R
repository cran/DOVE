#' Print the Primary Results of a dove() or dove2() Analysis
#'
#' Print the primary results of a dove() or dove2() analysis.
#'
#' @param x An DOVE object. The value object returned by a call to dove() or
#'   dove2()
#'
#' @param ... ignored
#'
#' @name print
#' @method print DOVE
#'
#' @returns No return value, called to display key results.
#'
#' @examples
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
#' print(x = result)
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
#' print(x = result)
#' 
#' @export
print.DOVE <- function(x, ...) {
  
  # 8/12/21: knots, gamma, covgamma, tau do not exist if type = 2

  # if (!is.null(attr(x = x, which = "knots"))) {
  #   attr(x = x, which = c("knots", "gamma", "covgamma", "tau", "type")) <- NULL
  # } else {
  #   attr(x = x, which = c("tau", "dSum", "type")) <- NULL
  # }
  
  if (!is.null(attr(x = x, which = "dSum"))) {
    attr(x = x, which = c("tau", "dSum", "type")) <- NULL
  } else {
    attr(x = x, which = "type") <- NULL
  }

  x <- unclass(x = x)
  print(x = x, ...)
}

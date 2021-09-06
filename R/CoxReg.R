# Internal function
#
# @param dt data.frame of time variables
# @param X matrix or NULL, covariates
# @param knots list of change point(s) of log hazard ratio
# @param constantVE TRUE: constant VE; FALSE: waning VE
# @param timePts numeric upper bound of intervals
# @param plots T/F, whether to plot curves of VE_a and VE_h or not

#' @include postProcess.R 
#' @importFrom stats pnorm
#' @import Rcpp
#' @importFrom stats sd
# 6/30/2021 extracted from iDOVE v1.2 with minor modification to
# conform to naming conventions of package and removal of right/left
CoxReg <- function(dt,  
                   X,  
                   knots,
                   constantVE,  
                   timePts,  
                   plots,
                   tau) {

  # number of participants
  n <- nrow(x = dt)

  if (is.null(x = X)) {
    # if no covariates provided define dummy variables
    p <- 0L
    varname <- NULL
    SD <- numeric(length = 0L)
    X <- matrix(data = 0.0, nrow = n, ncol = 1L)
  } else {
    p <- ncol(x = X)
    varname <- colnames(x = X)

    # standardize covariates X
    SD <- apply(X = X, MARGIN = 2L, FUN = sd, na.rm = TRUE)
    X <- scale(x = X, center = FALSE, scale = SD)
  }
  
  if (constantVE) {
    npc <- length(x = knots[[1L]])
  } else {
    npc <- length(x = knots[[1L]]) + 1L
  }

  eventTime <- dt$eventTime
  delta <- dt$eventStatus
  
  # sort the data by eventTime
  sorted_index <- order(eventTime)
  eventTime <- eventTime[sorted_index]
  delta <- delta[sorted_index]
  dt <- dt[sorted_index,,drop=FALSE]
  X <- X[sorted_index,,drop=FALSE]
    
  args <- list("beta"= rep(x = 0.0, times = p),
               "gamma" = rep(x = 0.0, times = npc), 
               "time" = eventTime,
               "delta" = delta,
               "X" = X,
               "W" = dt$entryTime,
               "S" = dt$vacTime, 
               "constantVE" = constantVE, 
               "threshold" = 0.0001,  
               "maxit" = 500)

  if (p == 0L) {
    funcCox <- "Cox_noX"
    Coxargs <- c("gamma", "time", "delta", "W", "S", 
                 "knots", "constantVE", "threshold", "maxit")
  } else {
    funcCox <- "Cox"
    Coxargs <- c("beta", "gamma", "time", "delta", "X", 
                 "W", "S", "knots", "constantVE", "threshold", 
                 "maxit")
  }

  # fit the standard Cox model

  if (length(x = knots) == 1L) {

    args[[ "knots" ]] <- knots[[ 1L ]]

    fit.cox <- tryCatch(expr = do.call(what = funcCox, args = args[ Coxargs ]),
                        error = function(e) {
                                  message(e$message)
                                  return( NULL )
                                })

    if (is.null(x = fit.cox)) stop("calculation aborted", call. = FALSE)

    ll0 <- fit.cox[[ 1L ]]
    covar <- fit.cox[[ 2L ]]
    theta <- fit.cox[[ 3L ]]

    finalknot <- 1L

    if (!anyNA(x = theta)) {
      message("Number of subjects: ", n) 
      message("log partial-likelihood at final estimates: ", 
              ifelse(test = abs(ll0) > 1.0, 
                     yes = round(x = ll0, digits = 2L),
                     no = round(x = ll0, digits = 4L)))
    } else {
      stop("NA values were produced in the standard Cox regression", 
           call. = FALSE)
    }

  } else {

    curll <- -Inf
    finalknot <- NA

    for (i in 1L:length(x = knots)) {


      args[[ "beta" ]] <- rep(x = 0.0, times = p)
      args[[ "gamma" ]] <- rep(x = 0.0, times = npc)

      args[[ "knots" ]] <- knots[[ i ]]

      fit.cox <- tryCatch(expr = do.call(what = funcCox, args = args[ Coxargs ]),
                          error = function(e) {
                                    message(e$message)
                                    return( NULL )
                                  })

      if (is.null(x = fit.cox)) next

      ll <- fit.cox[[ 1L ]]
      temptheta <- fit.cox[[ 3L ]]
      
      if (!anyNA(x = temptheta) & ll>curll) {
        finalknot <- i
        curll <- ll
        theta <- temptheta
        covar <- fit.cox[[ 2L ]]
      }
    }
    
    if (!is.na(x = finalknot)) {
      message("Day ", knots[[ finalknot ]], " (week ", knots[[ finalknot ]]/7, 
              ") was selected as the change point by AIC")
      message("Number of subjects: ", n) 
      message("log partial-likelihood at final estimates: ",
              ifelse(test = abs(curll) > 1.0, 
                     yes = round(x = curll, digits = 2L),
                     no = round(x = curll, digits = 4L)))
      
    } else {
      stop("all candidate change points produced NA values. ",
           "No change point was selected by AIC", call. = FALSE)
    }
  }

  knots <- knots[[ finalknot ]]

  res <- .postProcess(theta = theta, 
                      covMat = covar, 
                      plots = plots, 
                      SD = SD, 
                      varname = varname, 
                      constantVE = constantVE, 
                      tau = tau,
                      knots = knots,
                      timePts = timePts)

  res$changePts <- knots
  
  return( res )
}

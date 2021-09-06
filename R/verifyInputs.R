# input verification steps moved to new location to allow for use with dove2()
# added stopping condition for no events for vaccinated participants
# 6/30/2021 introduced for DOVE v1.7
#' @importFrom survival Surv is.Surv
#' @importFrom stats update.formula
.verifyInputs <- function(formula, data, timePts){

  data <- .verifyData(formula = formula, data = data)

  timePts <- .verifyTimePts(timePts = timePts)

  return( list("X" = data$X, 
               "data" = data$data, 
               "tau" = data$tau, 
               "timePts" = timePts) )

}


#' @importFrom stats terms update.formula 
#' @importFrom stats model.frame model.matrix model.response
#' @import survival 
#
# Note for developer: missing data are removed. X is either a matrix of
# the covariates or NULL indicating that there are no covariates. data
# is a data.frame containing {"entryTime", "eventTime", "eventStatus",
# "vacTime", "vacStatus"} Any NA vaccination times have been reset to tau
# tau is the ceiling of the maximum event time (considering only
# actual events time (event status = 1)).
.verifyData <- function(formula, data) {

  # separate vaccine object from formula
  # note for developer: extractVaccine() introduces an intercept for
  # Surv() ~ X component of formula
  formula <- .extractVaccine(formula = formula)
  vaccineForm <- formula$vaccineForm
  formula <- formula$formula

  # reset options to allow for keeping NA values 
  opt <- options()
  options(na.action = 'na.pass')
  on.exit(options(opt))
  
  # try to obtain the model.frame for Surv() ~ X component
  mf <- tryCatch(expr = stats::model.frame(formula = formula, data = data),
                 error = function(e) {
                   message("unable to obtain model.frame")
                   stop(e$message, call. = FALSE)
                 })
  
  # extract covariates (X is a matrix)
  # suppressing messages because already generated in creating mf
  X <- suppressMessages(expr = stats::model.matrix(object = mf, data = data))

  # remove intercept
  int <- attr(x = X, which = "assign") != 0L
  X <- X[,int, drop = FALSE]
  rownames(x = X) <- NULL

  if (ncol(x = X) == 0L) X <- NULL

  # dts will be a Surv object with columns "time" and "status"
  # again, messages suppressed because they will have already been
  # generated in the model.frame() call
  dts <- suppressMessages(expr = stats::model.response(data = mf))

  # if it is not a Surv object throw error
  if (!survival::is.Surv(x = dts)) {
    stop("the LHS of formula must be a Surv object", call. = FALSE)
  }

  # if it in interval Surv object throw error
  if (ncol(x = dts) != 2L) {
    stop("the Surv object must include event time and event status",  
         call. = FALSE)
  }

  dfs <- data.frame("eventTime" = dts[,1L], 
                    "eventStatus" = as.integer(x = round(x = dts[,2L], 
                                                         digits = 0L)))

  # try to obtain the model.frame for the vaccine() expression
  mfv <- tryCatch(expr = stats::model.frame(formula = vaccineForm, data = data),
                  error = function(e) {
                    message("unable to obtain model.frame")
                    stop(e$message, call. = FALSE)
                  })
  
  # dtv will be a vaccine() object with columns "entryTime", 
  # "vacTime", and "vacStatus"
  # again, messages suppressed because they will have already been
  # generated in the model.frame() call
  dtv <- suppressMessages(expr = stats::model.response(data = mfv))

  if (!is(object = dtv, class2 = "vaccine")) {
    stop("the RHS of formula must contain a vaccine() object", call. = FALSE)
  }

  dfv <- data.frame("entryTime" = dtv[,1L],
                    "vacTime" = dtv[,2L],
                    "vacStatus" = as.integer(x = round(x = dtv[,3L], 
                                                       digits = 0L)))

  df <- cbind(dfs, dfv)

  # remove any cases that have NA
  # NA may be in entry time indicating that time relationships were not
  # satisfied, or in event time and event status based on Surv()
  use <- stats::complete.cases(X, df)

  if (sum(use) == 0L) {
    stop("input checks result in all NA -- verify inputs", call. = FALSE)
  }

  if (sum(!use) > 0L) {
    df <- df[use,]
    if (!is.null(x = X)) X <- X[use,, drop = FALSE]
    message(sum(!use), " cases removed from the analysis due to NA values")
  }

  # ensure that event times >= entry times and that
  # event times >= vaccination times
  # vaccine() already tested to ensure that vacTime >= entryTime
  notOK <- {df$eventTime < {df$entryTime - 1e-8}} |
           {{df$eventTime < {df$vacTime - 1e-8}} & 
            {df$vacStatus == 1L} &
            {df$eventStatus == 1L}}

  violate <- sum(notOK)
  if (violate > 0L) {
    message(violate, 
            ifelse(test = violate > 1L, 
                   yes = " cases do not ", 
                   no = " case does not "),
            "satisfy required entry_time <= event_time and/or ",
            "vaccination_time <= event_time relationships; ",
            ifelse(test = violate > 1L, 
                   yes = "cases removed", 
                   no = "case removed"))
    df <-df[!notOK,]
    if (!is.null(x = X)) X <- X[!notOK,, drop = FALSE]
  } 

  if (sum({df$vacStatus == 1L} & {df$eventStatus == 1L}) == 0L) {
    stop("no vaccinated participants experienced an event; ", 
         "vaccine efficacy cannot be estimated", call. = FALSE)
  }

  # set tau to ceiling of maximum event time
  tau <- ceiling(x = max(df$eventTime[df$eventStatus == 1L]))

  message("tau = ", tau)
  if (tau <= 28.0) {
    warning("tau < 28 days (4 wks); ",
            "note: method expects time data to be on scale of days")
  }

  # reset Inf vaccination times
  df$vacTime[df$vacStatus == 0L] <- tau

  rownames(x = df) <- NULL

  return( list("X" = X, "data" = df, "tau" = tau) )
}


#' @importFrom stats terms reformulate
.extractVaccine <- function(formula) {

  cf <- as.character(x = formula)

  if (length(x = cf) != 3L) {
    stop("formula must include both LHS and RHS", call. = FALSE)
  }

  terms <- stats::terms(x = formula)

  termLabels <- attr(x = terms, which = "term.labels")

  vaccineObj <- which(x = grepl(pattern = 'vaccine(', 
                                x = termLabels, 
                                fixed = TRUE))

  if (length(x = vaccineObj) != 1L) {
    stop("unable to identify vaccine() object in formula", call. = FALSE)
  }

  response <- cf[2L]
  cov <- termLabels[-vaccineObj]
  vaccineObj <- termLabels[vaccineObj]

  if (length(x = cov) > 0L) {

    formula <- stats::reformulate(termlabels = cov, 
                                  response = response,
                                  intercept = TRUE)
  } else {
    formula <- stats::reformulate(termlabels = "1", 
                                  response = response,
                                  intercept = TRUE)

  }

  vaccineForm <- stats::reformulate(termlabels = "1", 
                                    response = vaccineObj,
                                    intercept = TRUE)

  return( list("formula" = formula, "vaccineForm" = vaccineForm) )  

}


.verifyTimePts <- function(timePts) {

  if (is.null(x = timePts)) {
    return( NULL )
  }

  if (!is.numeric(x = timePts) || !is.vector(x = timePts)) {
    stop("timePts must be a numeric vector or NULL", call. = FALSE)
  }
    
  if (length(x = timePts) == 0L) {
    return( NULL )
  } else if (any(timePts < 0.0)) {
    stop("all timePts must be positive", call. = FALSE)
  }

  return( timePts )
}



#' Durability of Vaccine Efficacy
#'
#' Estimates the potentially waning long-term efficacy of vaccines in
#'  randomized, placebo-controlled clinical trials with staggered
#'  enrollment of participants and sequential crossover of placebo recipients
#'
#' The information required for an analysis is 
#'   \describe{
#'     \item{Entry Time:}{Calendar time when the participant enters the 
#'       trial.}
#'     \item{Event Time:}{Calendar time when the participant experiences
#'       the clinical event of interest (e.g., symptomatic COVID-19) or their
#'       follow-up ends, whichever occurs first.}
#'     \item{Event Status:}{Binary indicator taking value 1 if the clinical
#'       event of interest occurs before the end of follow-up and 0 otherwise.}
#'     \item{Vaccination Status:}{Binary indicator taking value 1 if
#'       vaccination occurs before the end of follow-up and 0 otherwise.}
#'     \item{Vaccination Time:}{Calendar time when vaccination takes place,
#'       with an arbitrary value if the participant is not vaccinated.}
#'     \item{Covariates:}{Baseline covariates (e.g., priority group, age, 
#'       ethnicity).}
#'    }
#'
#' Note that all times are to be specified relative to the start of the trial.
#' Thus, for individuals that received vaccination, 
#' entry_time <= vaccination_time <= event_time. And for
#' individuals that did not receive vaccination, 
#' entry_time <= event_time; for these participants,
#' vaccination_time can take any value.
#'
#' In plots, it is assumed that times are specified in units of days, though
#' this is not verified nor required.
#' 
#' The general structure of the formula input is
#'   \preformatted{
#'   Surv(event_time, event_status) ~ covariates + 
#'     vaccine(entry_time, vaccination_time, vaccination_status)
#'   }
#'
#' The response variable must be a survival object as returned by the
#' 'Surv()' function of package \pkg{survival}, where event_time is the
#' follow up time (formal argument 'time') and event_status is the status
#' indicator input (formal argument 'event'). Specifically, 
#' \preformatted{Surv(time = event_time, event = event_status)}
#'
#' The vaccination and entry_time information must be specified through function 
#' 'vaccine()'. Specifically, 
#' \preformatted{vaccine(entry_time, vaccination_status, vaccination_time)}
#' For participants that did not receive the vaccine, vaccination_time
#' can take any value. For individuals
#' that received vaccination, if 
#' vaccination_time > event_time, or vaccination_time < entry_time, the 
#' case will be removed from the analysis and a message will be generated.
#'
#' @rdname dove
#' @name dove
#' 
#' @references Lin, DY, Zeng, D, and Gilbert, PB (2021). Evaluating the 
#'   long-term efficacy of COVID-19 vaccines. 
#'   doi: https://doi.org/10.1101/2021.01.13.21249779.
#'
#' @param formula A formula object, with the response on the left hand side of a
#'   '~' operator, and the covariates and vaccine() function on the right.  
#'   The response must be a survival object as returned by the 'Surv'
#'   function of the \pkg{survival} package. See Details for further information.
#'   The vaccine() function must be used to specify the entry time and
#'   vaccination information. See ?vaccine and Details for further
#'   information. 
#' 
#' @param data A data.frame object. The data.frame in which to interpret the
#'   variable names in formula. Must contain the entry time, the event time,
#'   the event status, the vaccination status,
#'   the vaccination time, and the covariates. See Details.
#'   
#' @param plots A logical object. If TRUE, plots of the estimated
#'   curve of vaccine efficacy in reducing attack rate and
#'   the estimated curve of vaccine efficacy in reducing the hazard rate
#'   will be automatically generated. If FALSE, plots will not be
#'   generated, but the data are available through the returned value object.
#'
#' @param timePts A numeric vector object or NULL. The endpoints of the time 
#'   periods for which the vaccine efficacy in reducing the attack rate is
#'   to be shown. If NULL, a default sequence of 60 day intervals is used.
#'
#' @param bandwidth A numeric vector object. A tuning parameter for the 
#'    bandwidth used for kernel estimation of the vaccine efficacy in 
#'    reducing the hazard rate; this input is ignored if plots=FALSE.
#'
#' @returns A list containing
#'   \item{covariates}{A matrix containing the estimated hazard ratio of each
#'     covariate, together with the (estimated) standard error, the 95\%
#'     confidence intervals, and the two-sided p-value for testing no covariates
#'     effect.}
#'   \item{vaccine}{A list containing two elements. The first element is
#'     the matrix containing the estimates of the vaccine efficacy in
#'     reducing the attack rate at all observed event 
#'     times (VE_a), together with the 95\% 
#'     confidence intervals, as well as the vaccine efficacy in
#'     reducing the hazard rate at these times (VE_t). These results 
#'     will be shown in graphical form if input plots = TRUE.
#'     The second element is the matrix containing the estimates of vaccine
#'     efficacy in reducing the attack rate over successive time periods,
#'     together with the 95\% confidence intervals.
#'     }
#' 
#' @export
#' @import methods
#' 
#' @include vaccine.R vcFit.R
#' @importFrom graphics axis legend lines plot plot.new
#' @importFrom stats model.response update.formula complete.cases 
#' @importFrom stats stepfun plot.stepfun
#' @importFrom grDevices dev.new
#' @importFrom survival Surv is.Surv
#'
#' @examples
#'
#' data(doveData)
#'
#' set.seed(1234)
#'
#' ind <- sample(1:nrow(x = doveData), 2500, FALSE)
#'
#' # NOTE: This sample size is chosen for example only -- larger data sets
#' # should be used.
#' # See the vignette for a full analysis of the doveData dataset
#'
#' dove(formula = Surv(event.time, event.status) ~ priority + sex + 
#'                 vaccine(entry.time, vaccine.status, vaccine.time),
#'      data = doveData[ind,])

dove <- function(formula, data, plots = TRUE, timePts = NULL, bandwidth = NULL) {

  if (missing(x = formula)) {
    stop("a formula argument must be provided", call. = FALSE)
  }

  if (missing(x = data)) {
    stop("a data argument must be provided", call. = FALSE)
  }

  # reset options to allow for keeping na values
  opt <- options()
  options(na.action = 'na.pass')
  on.exit(options(opt))

  # add intercept from model if not provided
  # this was added 6/2/21 to ensure that factors are handled properly
  if (attr(x = stats::terms(x = formula), which = "intercept") == 0L) {
    formula = update.formula(old = formula, new = .~. +1)
  }

  # try to obtain the model.frame
  mf <- tryCatch(expr = stats::model.frame(formula = formula, data = data),
                 error = function(e) {
                           message("unable to obtain model.frame")
                           stop(e$message, call. = FALSE)
                          })

  # extract covariates
  X <- suppressMessages(stats::model.matrix(object = mf, data = data))
  # 6/2/21 remove intercept
  int <- attr(x = X, which = "assign") != 0L
  X <- X[,int, drop = FALSE]

  # identify the columns that correspond to the returns returns by 
  # vaccine()

  lbl <- attr(x = stats::terms(x = formula), which = "term.labels")

  lbl1 <- paste0(lbl, "vaccination_time")
  lbl2 <- paste0(lbl, "vaccination_status")
  lbl3 <- paste0(lbl, "entry_time")

  if (!any(lbl1 %in% colnames(x = X)) || 
      !any(lbl2 %in% colnames(x = X)) || 
      !any(lbl3 %in% colnames(x = X))) {
    stop("the RHS of formula did not contain an appropriate vaccine() object",
         call. = FALSE)
  }

  i1 <- colnames(x = X) %in% lbl1
  i2 <- colnames(x = X) %in% lbl2
  i3 <- colnames(x = X) %in% lbl3

  vacTime <- X[,i1]
  vacStatus <- X[,i2]
  entryTime <- X[,i3]
  # 6/2/21 added drop=FALSE to properly handle case when only 1 covariate
  # is in the model.
  X <- X[,-c(which(i1),which(i2),which(i3)),drop = FALSE]

  if (ncol(x = X) == 0L) X <- NULL
  
  # dt will be a Surv object with columns "time" and "status"
  dt <- suppressMessages(stats::model.response(data = mf))

  if (!is.Surv(x = dt)) {
    stop("the LHS of formula must be a Surv object", call. = FALSE)
  }

  if (ncol(x = dt) != 2L) {
    stop("the Surv object must include event time, and", 
         " event status", call. = FALSE)
  }

  eventTime <- dt[,1L]
  eventStatus <- dt[,2L]

  # remove any cases that have NA in the 
  # entry_time, event_time, event_status, covariates, or vacStatus
  use <- stats::complete.cases(cbind(X,entryTime,eventTime,eventStatus,vacStatus))

  if (all(!use, na.rm = TRUE)) {
    stop("input checks result in all NA -- verify inputs",
         call. = FALSE)
  }

  if (any(!use, na.rm = TRUE)) {
    entryTime <- entryTime[use]
    eventTime <- eventTime[use]
    eventStatus <- eventStatus[use]
    if (!is.null(x = X)) X <- X[use,, drop = FALSE]
    vacTime <- vacTime[use]
    vacStatus <- vacStatus[use]
  } 

  # ensure that event times are after entry times and that
  # vaccination times are before event times
  # vaccine() already tested to ensure that vacTime >= entryTime
  tst <- {{eventTime < entryTime} | {vacTime > eventTime}}

  if (any(tst, na.rm = TRUE)) {
    entryTime <- entryTime[!tst]
    eventTime <- eventTime[!tst]
    eventStatus <- eventStatus[!tst]
    if (!is.null(x = X)) X <- X[!tst,, drop = FALSE]
    vacTime <- vacTime[!tst]
    vacStatus <- vacStatus[!tst]
  } 

  if (sum(!use, na.rm = TRUE) > 0L || sum(tst, na.rm = TRUE) > 0L) {
    message(sum(!use | tst, na.rm = TRUE), 
            " cases removed from the analysis due to NA values")
  }

  # set tau to epsilon above the maximum event time
  tau <- max(eventTime[eventStatus == 1L], na.rm = TRUE) + 1e-8
  message("tau set to ", round(x = tau, digits = 4L))

  # remove NA vaccination times for convenience
  vacTime[vacStatus == 0L] <- tau

  # verify timePts
  if (!is.null(x = timePts)) {
    if (!is.numeric(x = timePts)) {
      stop("timePts must be a numeric vector or NULL", call. = FALSE)
    }

    if (length(x = timePts) <= 0L) {
      message("timePts has zero length; default values will be used")
      timePts <- NULL
    } else if (any(timePts < 0.0)) {
        stop("all timePts must be positive", call. = FALSE)
    }

  } 

  if (is.null(x = timePts)) {
    timePts <- c(seq(from = 60L, to = tau, by = 60L), tau)
    if (any(diff(x = timePts) < 60L)) {
      nt <- length(x = timePts)
      timePts <- timePts[-nt]
    }
  }

  timePts <- sort(x = unique(x = timePts))

  if (is.null(x = bandwidth)) bandwidth = 0.5

  if (!is.numeric(x = bandwidth)) {
    stop("bandwidth must be numeric", call. = FALSE)
  }

  if (length(x = bandwidth) > 1L) {
    message("using only first bandwidth")
    bandwidth <- bandwidth[1L]
  }

  res <- .vcFit(R = entryTime, 
                Y = eventTime,  
                Delta = as.integer(x = round(x = eventStatus, digits = 0L)),  
                D = as.integer(x = round(x = vacStatus, digits = 0L)),  
                S = vacTime,  
                X = X,  
                tau = tau,
                timePts = timePts)

  class(res) <- "doveObj"

  attr(res, "tau") <- tau

  Delta <- as.integer(x = round(x = eventStatus, digits = 0L))
  D <- as.integer(x = round(x = vacStatus, digits = 0L))
 
  attr(res, "dSum") <- sum(Delta[D == 1L] == 1L)
  
  if (plots) plot(x = res, bandwidth = bandwidth)

  # until the plot feature is exported -- no reason to class object
  attr(x = res, "tau") <- NULL
  attr(x = res, "dSum") <- NULL
  res <- unclass(x = res)


  return( res )
}

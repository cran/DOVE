#' Durability of Vaccine Efficacy
#'
#' Estimates the potentially waning long-term efficacy of vaccines in
#'  randomized, placebo-controlled clinical trials with staggered
#'  enrollment of participants and sequential crossover of placebo recipients.
#'  The hazard ratio for the vaccine effect is a nonparametric 
#'  function of time.
#'
#' In dove(), the hazard ratio for the vaccine effect is a nonparametric 
#'  function of time for which
#'  confidence intervals are not provided. The sister function, dove2(), 
#'  assumes that the log hazard ratio for the vaccine effect is a piecewise 
#'  linear function of time; 
#'  it provides more precise and more stable estimates of vaccine efficacy on 
#'  the hazard rate and includes proper confidence intervals.
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
#'       with an arbitrary non-negative value if the participant is not vaccinated.}
#'     \item{Covariates:}{Baseline covariates (e.g., priority group, age, 
#'       ethnicity).}
#'    }
#'
#' Note that all times are to be specified relative to the start of the trial
#' and are specified in units of days. Thus, for individuals that received 
#' vaccination, entry_time <= vaccination_time <= event_time. And for
#' individuals that did not receive vaccination, 
#' entry_time <= event_time; for these participants,
#' vaccination_time can take any non-negative value (including NA or Inf).
#'
#' The general structure of the formula input is
#'   \preformatted{
#'   Surv(event_time, event_status) ~ covariates + 
#'     vaccine(entry_time, vaccination_time, vaccination_status)
#'   }
#'
#' The response variable must be a survival analysis object as returned by the
#' 'Surv()' function of package \pkg{survival}, where event_time is the
#' observation time (formal argument 'time') and event_status is the status
#' indicator input (formal argument 'event'). Specifically, 
#' \preformatted{Surv(time = event_time, event = event_status)}
#'
#' The covariates can include categorical variables, for which
#' all other categories are compared to the first category.
#'
#' The vaccination and entry_time information must be specified through function 
#' 'vaccine()'. Specifically, 
#' \preformatted{vaccine(entry_time, vaccination_status, vaccination_time)}
#' For participants that did not receive the vaccine, vaccination_time
#' can take any non-negative value (including NA or Inf). For individuals
#' that received vaccination, if 
#' vaccination_time > event_time, or vaccination_time < entry_time, the 
#' case will be removed from the analysis and a message will be generated.
#'
#' @rdname dove
#' @name dove
#' 
#' @references Lin, DY, Zeng, D, and Gilbert, PB (2021). Evaluating the 
#'   long-term efficacy of COVID-19 vaccines. Clinical Infectious Diseases, 
#'   ciab226, https://doi.org/10.1093/cid/ciab226.
#'
#' @param formula A formula object, with the response on the left-hand side of a
#'   '~' operator, and the covariates and vaccine() function on the right-hand side.  
#'   The response must be a survival object as returned by the 'Surv'
#'   function of the \pkg{survival} package. See Details for further information.
#'   The vaccine() function must be used to specify the entry time and
#'   vaccination information. See ?vaccine and Details for further
#'   information. 
#' 
#' @param data A data.frame object. The data.frame in which to interpret the
#'   variable names in formula. Must contain the entry time, the event time,
#'   the event status, the vaccination status,
#'   the vaccination time, and any covariates. All time variables must
#'   be provided in units of days. See Details.
#'   
#' @param plots A logical object. If TRUE, plots of the estimated
#'   curve of vaccine efficacy in reducing attack rate and
#'   the estimated curve of vaccine efficacy in reducing the hazard rate
#'   will be automatically generated. If FALSE, plots will not be
#'   generated, but the data are available through the returned value object.
#'
#' @param timePts A numeric vector object or NULL. The endpoints of the time 
#'   periods for which the vaccine efficacy in reducing the attack rate is
#'   to be shown. If NULL, a default sequence of 60-day intervals is used.
#'   If tau < 60 days, this input must be provided.
#'
#' @param bandwidth A numeric object. Tuning parameter for the 
#'    bandwidth used for kernel estimation of the vaccine efficacy in 
#'    reducing the hazard rate; this input is ignored if plots=FALSE.
#'
#' @returns An S3 object of class DOVE containing a list with elements
#'   \item{covariates}{A matrix containing the estimated hazard ratio of each
#'     covariate, together with the (estimated) standard error, the 95\%
#'     confidence interval, and the two-sided p-value for testing no covariate
#'     effect.}
#'   \item{vaccine}{A list containing two elements. The first element is
#'     the matrix containing the estimates of the vaccine efficacy in
#'     reducing the attack rate at all observed event 
#'     times (VE_a), together with the 95\% 
#'     confidence intervals, as well as the vaccine efficacy in
#'     reducing the hazard rate at these times (VE_h). These results 
#'     will be shown in graphical form if input plots = TRUE.
#'     The second element is the matrix containing the estimates of vaccine
#'     efficacy in reducing the attack rate over successive time periods,
#'     together with the 95\% confidence intervals.
#'     }
#'
#' @export
#' @import methods
#' @importFrom RcppArmadillo armadillo_throttle_cores
#' 
#' @include verifyInputs.R vcFit.R
#'
#' @examples
#' \dontshow{
#'   RcppArmadillo::armadillo_throttle_cores(2)
#' }
#' data(doveData)
#'
#' set.seed(1234)
#'
#' ind <- sample(1:nrow(x = doveData), 500, FALSE)
#'
#' # NOTE: This sample size is chosen for example only -- larger data sets
#' # should be used.
#' # See the vignette for a full analysis of the doveData dataset
#'
#' dove(formula = Surv(event.time, event.status) ~ priority + sex + 
#'                 vaccine(entry.time, vaccine.status, vaccine.time),
#'      data = doveData[ind,])
# 6/30/2021 modified to shift input verification to external function
# used by both dove and dove2. S3 class set as DOVE
dove <- function(formula, 
                 data, 
                 plots = TRUE, 
                 timePts = NULL, 
                 bandwidth = NULL) {

  if (missing(x = formula)) {
    stop("a formula argument must be provided", call. = FALSE)
  }

  if (missing(x = data)) {
    stop("a data argument must be provided", call. = FALSE)
  }

  # matched call to include with returned object
  cl <- match.call()
  
  # process inputs to obtain time and covariate information and to
  # verify timePts
  inputs <- .verifyInputs(formula = formula, data = data, timePts = timePts)

  # notify user if default timePts will be used and set timePts
  if (is.null(x = inputs$timePts)) {
    message("timePts not given; default values will be used")

    if (inputs$tau > 60L) {

      timePts <- c(seq(from = 60L, to = inputs$tau, by = 60L), inputs$tau)

      if (any(diff(x = timePts) < 60L)) {
        nt <- length(x = timePts)
        timePts <- timePts[-nt]
      }
    } else {
      stop("tau < 60 days; define timePts on input", call. = FALSE)
    }
  } else {
    timePts <- inputs$timePts
  }

  timePts <- sort(x = unique(x = timePts))

  if (is.null(x = bandwidth)) bandwidth = 0.3

  if (!is.numeric(x = bandwidth)) {
    stop("bandwidth must be numeric", call. = FALSE)
  }

  res <- .vcFit(R = inputs$data$entryTime, 
                Y = inputs$data$eventTime,  
                Delta = inputs$data$eventStatus,  
                D = inputs$data$vacStatus,  
                S = inputs$data$vacTime,  
                X = inputs$X,  
                tau = inputs$tau,
                timePts = timePts)

  res[["call"]] <- cl

  class(x = res) <- "DOVE"

  attr(x = res, which = "tau") <- inputs$tau

  Delta <- inputs$data$eventStatus
  D <- inputs$data$vacStatus
 
  attr(x = res, which = "dSum") <- sum(Delta[D == 1L] == 1L)
  attr(x = res, which = "type") <- 1L
  
  if (plots) plot(x = res, bandwidth = bandwidth)

  return( res )
}

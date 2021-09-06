# Internal function
#
# @param R numeric vector, entry times
# @param Y numeric vector, event times
# @param Delta integer vector, event status
# @param D integer vector, vaccination status
# @param S numeric vector, vaccination times
# @param X matrix or NULL, covariates
# @param tau numeric
# @param timePts numeric upper bound of intervals
#' @importFrom stats pnorm quantile
# 7/30/2021 -- constrain VE_a and VE_h to be 0 at day 0 (lines 394-403) 
#' @include loc.R scores.R postProcessDove.R
.vcFit <- function(R, Y, Delta, D, S, X, tau, timePts) {

  # number of basis function
  m <- 20L

  # number of participants in study
  n <- length(x = R)

  # number of covariates (if any)
  if (is.null(x = X)) {
    np <- 0L
  } else {
    np <- ncol(x = X)
  }

  # times are taken as the 1/m, 2/m, ..., m-1/m disease times
  # type = 5 appears to match matlab percentil method
  # {m-1}
  times <- stats::quantile(x = Y[Delta == 1L], probs = {1L:{m-1L}}/m, type = 5L)

  # {m}
  timesLower <- unname(obj = c(0, times))

  # {m}
  timesUpper <- unname(obj = c(times, tau))

  # components of score1 that need only be calculated once
  score1Prep <- .score1Prep(X = X,
                            Delta = Delta,
                            D = D,
                            YStar = pmin(Y, S),
                            timesLower = timesLower,
                            timesUpper = timesUpper,
                            R = R)

  # vaccinated
  vac <- D == 1L

  # {nVac}
  diseaseTimeAfterVac <- Y[vac] - S[vac]

  # events for vaccinated participants
  # {nVac}
  DeltaVac <- Delta[vac]

  # ensure that the last time step does not reflect an event
  DeltaVac[diseaseTimeAfterVac > max(diseaseTimeAfterVac) - 1e-8] <- 0L 

  # unique event times for vaccinated participants
  eventTimesAfterVac <- diseaseTimeAfterVac[DeltaVac == 1L]

  # indicator of time since vaccination >= event times lower bounds
  # Y1-S1 >= T1 Y1-S1 >= T2 ... Y1-S1 >= Tx
  # Y2-S2 >= T1 Y2-S2 >= T2 ... Y2-S2 >= Tx
  # ...
  # Yn-Sn >= T1 Yn-Sn >= T2 ... Yn-Sn >= Tx
  # {nVac x nEventsVac}
  indVac <- outer(X = diseaseTimeAfterVac, 
                  Y = eventTimesAfterVac-1e-8, 
                  FUN = ">") 

  loc <- .loc(S = S[vac],
              eventTimes = eventTimesAfterVac,
              timesLower = timesLower,
              timesUpper = timesUpper)

  # add np to simplify pulling correct beta
  loc <- loc + np

  # components of score that need only be calculated once
  score2Prep <- .score2Prep(vac = vac,
                            Delta = Delta,
                            X = X,
                            Y = Y,
                            timesLower = timesLower,
                            timesUpper = timesUpper)

  if (is.null(x = X)) {
    itr <- tryCatch(expr = iteration_NoX(score1Prep$dt, score1Prep$DZBBy,
                                         indVac, loc, score2Prep, m),
                      error = function(e) {
                                stop(e$message, call. = FALSE)
                               })
  } else {
    itr <- tryCatch(expr = iteration(score1Prep$dt, score1Prep$DZBBy,
                                     X, indVac, loc, score2Prep,
                                     m, vac*1L),
                      error = function(e) {
                                stop(e$message, call. = FALSE)
                               })
  }

  success <- itr[[ 1L ]]
  beta <- itr[[ 2L ]]
  s1 <- list("score" = itr[[ 3L ]], "dscore" = itr[[ 4L ]])
  s2 <- list("score" = itr[[ 5L ]], "dscore" = itr[[ 6L ]])

  message("method converged after ", success, " iterations")

  outputVaccine <- .postProcessDove(X = X,
                                    beta = beta,
                                    loc = loc,
                                    DeltaVac = DeltaVac,
                                    s1 = s1,
                                    s2 = s2,
                                    diseaseTimeAfterVac = diseaseTimeAfterVac,
                                    eventTimesAfterVac = eventTimesAfterVac,
                                    vac = vac,
                                    timePts = timePts,
                                    indVac = indVac)

  return( outputVaccine )
}

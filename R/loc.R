# note to developer: anything within 1e-8 is deemed "equal"
.loc <- function(S, eventTimes, timesLower, timesUpper) {

  nVac <- length(x = S)
  nEventTimes <- length(x = eventTimes)
  m <- length(x = timesLower)

  # all possible combinations of time to vaccination and time to event
  # after vaccination
  # {nVac x nEventsVac}
  Y2S <- outer(X = S, Y = eventTimes, FUN = "+")

  # identify the eventTimesVac bin in which Y2S falls
  loc <- matrix(data = 0L, nrow = nVac, ncol = nEventTimes)

  for (k in 1L:m) {
    loc <- loc + k*{Y2S > {timesLower[k]+1e-8}} * {Y2S < {timesUpper[k]+1e-8}}
  }

  return( loc )
}


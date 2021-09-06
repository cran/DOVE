.score1Prep <- function(X, Delta, D, YStar, timesLower, timesUpper, R) {

  n <- length(x = Delta)
  m <- length(x = timesUpper)

  # {m}
  dTimes <- timesUpper - timesLower

  # {n x m}
  mdTimes <- matrix(data = dTimes, nrow = n, ncol = m, byrow = TRUE)

  # rows correspond to YStar; columns to lower bound of each time step
  # matrix is an indicator matrix
  # Y1 > L1  Y1 > L2  Y1 > L3 ... Y1 > Lm
  # Y2 > L1  Y2 > L2  Y2 > L3 ... Y2 > Lm
  # ...
  # Yn > L1  Yn > L2  Yn > L3 ... Yn > Lm
  # {n x m}
  YStarLower <- outer(X = YStar, Y = {timesLower+1e-8}, FUN = ">")

  # rows correspond to YStar, columns to upper bound of each time step
  # matrix is an indicator matrix
  # Y1 <= U1  Y1 <= U2  Y1 <= U3 ... Y1 <= Um
  # Y2 <= U1  Y2 <= U2  Y2 <= U3 ... Y2 <= Um
  # ...
  # Yn <= U1  Yn <= U2  Yn <= U3 ... Yn <= Um
  # {n x m}
  YStarUpper <- outer(X = YStar, Y = {timesUpper+1e-8}, FUN = "<")

  # 1 = participant experienced disease and was not vaccinated
  # {n}
  diseaseNoVac <- Delta * (1L - D)

  # risk factors concatenated with matrix indicating the time bin in which the
  # response falls zeroed out for participants that were vaccinated and for
  # participants that were not vaccinated but are not diseased
  # {n x {np + m}}
  DZBBy <- diseaseNoVac*cbind(X, YStarLower*YStarUpper)

  # dt for integration from 0 to Y
  # {n x m}
  dtY <- YStarLower * pmin(outer(X = YStar, Y = timesLower, FUN = "-"), 
                           mdTimes)

  # dt for integration from 0 to T0
  # {n x m}
  RLower <- outer(X = R, Y = {timesLower-1e-8}, FUN = ">")
  dtEntryTime <- RLower*
                 pmin(outer(X = R, Y = timesLower, FUN = "-"), mdTimes)

  # dt for integration from T0 to Y
  # {n x m}
  dt <- dtY - dtEntryTime

  return( list("DZBBy" = unname(obj = DZBBy), "dt" = unname(obj = dt)) )
}


.score2Prep <- function(vac, Delta, X, Y, timesLower, timesUpper) {

  # risk factors concatenated with indicator matrix showing the time bin
  # in which the response falls
  # taken only for individuals that were vaccinated and experienced symptoms
  # {nEventsVac x {np+m}}
  ind <- vac & {Delta == 1L}

  ZBBy2 <- cbind(X[ind,,drop=FALSE], 
                 outer(X = Y[ind], Y = timesLower, FUN = ">")*
                 outer(X = Y[ind], Y = timesUpper+1e-8, FUN = "<"))

  return( ZBBy2 )

}

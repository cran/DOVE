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

  # YStar = (1-D)Y + DS where D = I(S <= Y)
  # {n}
  YStar <- pmin(Y, S)

  # 1 = participant experienced disease and was not vaccinated
  # {n}
  diseaseNoVac <- Delta * (1L - D)

  # times are taken as the 1/m, 2/m, ..., m-1/m disease times
  # type = 5 appears to match matlab percentil method
  # {m-1}
  times <- stats::quantile(x = Y[Delta == 1L], probs = {1L:{m-1L}}/m, type = 5L)

  # {m}
  timesLower <- c(0, times)

  # {m}
  timesUpper <- c(times, tau)

  # {m}
  dTimes <- timesUpper - timesLower

  # {n x m}
  mdTimes <- matrix(data = dTimes, nrow = n, ncol = m, byrow = TRUE)

  # row correspond to YStar; columns to lower bound of each time step
  # matrix is an indicator matrix
  # Y1 > L1  Y1 > L2  Y1 > L3 ... Y1 > Lm
  # Y2 > L1  Y2 > L2  Y2 > L3 ... Y2 > Lm
  # ...
  # Yn > L1  Yn > L2  Yn > L3 ... Yn > Lm
  # {n x m}
  YStarLower <- outer(X = YStar, Y = timesLower, FUN = ">")

  # row correspond to YStar, columns to upper bound of each time step
  # matrix is an indicator matrix
  # Y1 <= U1  Y1 <= U2  Y1 <= U3 ... Y1 <= Um
  # Y2 <= U1  Y2 <= U2  Y2 <= U3 ... Y2 <= Um
  # ...
  # Yn <= U1  Yn <= U2  Yn <= U3 ... Yn <= Um
  # {n x m}
  YStarUpper <- outer(X = YStar, Y = timesUpper, FUN = "<=")

  # dt for integration from 0 to Y
  # {n x m}
  dtY <- YStarLower * pmin(outer(X = YStar, Y = timesLower, FUN = "-"), 
                           mdTimes)

  # risk factors concatenated with matrix indicating the time bin in which the
  # response falls zeroed out for participants that were vaccinated and for
  # participants that were not vaccinated but are not diseased
  # {n x {np + m}}
  DZBBy <- diseaseNoVac*cbind(X, YStarLower*YStarUpper)

  # dt for integration from 0 to T0
  # {n x m}
  RLower <- outer(X = R, Y = timesLower, FUN = ">=")
  dtEntryTime <- RLower*
                 pmin(outer(X = R, Y = timesLower, FUN = "-"), mdTimes)

  # dt for integration from T0 to Y
  # {n x m}
  dt <- dtY - dtEntryTime

  # vaccinated
  vac <- D == 1L

  # totoal number of participants vaccinated
  nVac <- sum(vac)

  # {nVac}
  diseaseTimeAfterVac <- Y[vac] - S[vac]

  # events for vaccinated participants
  # {nVac}
  DeltaVac <- Delta[vac]

  # ensure that the last time step does not reflect an event
  DeltaVac[diseaseTimeAfterVac >= max(diseaseTimeAfterVac) - 1e-8] <- 0L 

  # unique event times for vaccinated participants
  eventTimesAfterVac <- diseaseTimeAfterVac[DeltaVac == 1L]

  # number of event times for vaccinated particpants
  nEventsVac <- length(x = eventTimesAfterVac)

  # indicator of time since vaccination >= event times lower bounds
  # Y1-S1 >= T1 Y1-S1 >= T2 ... Y1-S1 >= Tx
  # Y2-S2 >= T1 Y2-S2 >= T2 ... Y2-S2 >= Tx
  # ...
  # Yn-Sn >= T1 Yn-Sn >= T2 ... Yn-Sn >= Tx
  # {nVac x nEventsVac}
  indVac <- outer(X = diseaseTimeAfterVac, Y = eventTimesAfterVac, FUN = ">=") 

  # all possible combinations of time to vaccination and time to event
  # after vaccination
  # {nVac x nEventsVac}
  Y2S <- outer(X = S[vac], Y = eventTimesAfterVac, FUN = "+")

  # identify the evenTimesVac bin in which Y2S falls
  loc <- matrix(data = 0L, nVac, nEventsVac)

  for (k in 1L:m) {
    loc <- loc + k*{Y2S > timesLower[k]} * {Y2S <= timesUpper[k]}
  }
  # add np to simplify pulling correct beta
  loc <- loc + np

  # risk factors concatenated with indicator matrix showing the time bin
  # in which the response falls
  # taken only for individuals that were vaccinated and experienced symptoms
  # {nEventsVac x {np+m}}
  ind <- vac & Delta == 1L

  ZBBy2 <- cbind(X[ind,,drop=FALSE], 
                 outer(X = Y[ind], Y = timesLower, FUN = ">")*
                 outer(X = Y[ind], Y = timesUpper, FUN = "<="))

  oldbeta <- numeric(length = np+m)
  iter <- 0L 
  error <- 1.0 
  maxiter <- 500L
  epsilon <- 0.0001

  dscore1 <- matrix(data = 0.0, nrow = np+m, ncol = np+m)
  score1 <- matrix(data = 0.0, nrow = n, ncol = np+m)


  while ({error > epsilon} && {iter < maxiter}) {

    # {n}
    if (np > 0L) {
      zbeta <- drop(x = X %*% oldbeta[1L:np])
    } else {
      zbeta <- numeric(length = n)
    }

    # {n x m}
    EXbeta <- exp(x = outer(X = zbeta, 
                            Y = oldbeta[{np+1L}:{np+m}],  
                            FUN = "+"))

    # {n x m}
    exbt <- EXbeta*dt

    for (k in 1L:np) {
      # {n x m}
      tmp <- X[,k]*exbt

      score1[,k] <- DZBBy[,k] - rowSums(x = tmp)

      for (j in k:np) {
        dscore1[k,j] <- -sum(X[,j]*tmp)
        dscore1[j,k] <- dscore1[k,j]
      }
      for (j in {np+1L}:{np+m}) {
        dscore1[k,j] <- -sum(tmp[,j-np])
        dscore1[j,k] <- dscore1[k,j]
      }
    }

    for (k in {np+1L}:{np+m}) {
      score1[,k] <- DZBBy[,k] - exbt[,k-np]
    }

    diag(x = dscore1)[{np+1L}:{np+m}] <- -colSums(x = exbt)

    ####
    
    dscore2 <- matrix(data = 0.0, nrow = np+m, ncol = np+m)
    score2 <- matrix(data = 0.0, nrow = nEventsVac, ncol = np+m)

    # {nVac}
    zbeta2 <- zbeta[vac]

    # {nVac x nEventsVac}
    tmp <- matrix(data = oldbeta[loc], ncol = ncol(x = loc))
    tmp[loc<=np] <- 0L
    EXbeta <- exp(x = zbeta2 + tmp)

    # {nVac x nEventsVac}
    iexb <- indVac*EXbeta

    # {nEventsVac}
    denom <- colSums(x = iexb)

    # {np} each {nEventsVac}
    uuu <- list()
    k <- 1L
    while (k <= np) {
      uuu[[ k ]] <- colSums(iexb*X[vac,k])
      k <- k + 1L
    }

    # {nEventsVac x m}
    ttt <- sapply(X = {np+1L}:{np+m},
                  FUN = function(x,y,z){colSums(y*{z==x})},
                  y = iexb, z = loc)

    k <- 1L
    while (k <= np) {
      # d exp(Z beta + gamma) / d beta_k
      # {nVac x nEventsVac}
      tempk <- iexb*X[vac,k]

      # sum_{i=1}^n d exp(Z beta + gamma) / d beta_k
      # {nEventsVac}
      nomk <- uuu[[ k ]]

      score2[,k] <- ZBBy2[,k] - nomk/denom

      for (j in k:np) {
        # {nEventsVac}
        nomj <- uuu[[ j ]] 

        # {nVac x nEventsVac}
        tempjk <- tempk*X[vac,j]
        # {nEventsVac}
        nomjk <- colSums(x = tempjk)

        dscore2[k,j] <- -sum(nomjk/denom - nomj*nomk/{denom^2})
        dscore2[j,k] <- dscore2[k,j]
      }

      for (j in {np+1L}:{np+m}) {
        # {nEventsVac}
        nomj <- ttt[,j-np]

        # {nEventsVac}
        nomjk <- colSums(x = tempk*{loc==j})

        dscore2[k,j] <- -sum(nomjk/denom - nomj*nomk/{denom^2})
        dscore2[j,k] <- dscore2[k,j]
      }
      k <- k + 1L
    }

    for (k in {np+1L}:{np+m}) {
      # d exp(Z beta + gamma) / d gamma_l
      # {nVac x nEventsVac}
      tempk <- iexb*{loc==k}

      # sum_{i=1}^n d exp(Z beta + gamma) / d gamma_l
      # {nEventsVac}
      nomk <- ttt[,k-np]

      score2[,k] <- ZBBy2[,k] - nomk/denom

      dscore2[k,k] <- -sum(ttt[,k-np]/denom - ttt[,k-np]*ttt[,k-np]/{denom^2})

      j <- k + 1L
      while( j <= {np+m} ) {
        dscore2[k,j] <- sum(ttt[,j-np]*ttt[,k-np]/{denom^2})
        dscore2[j,k] <- dscore2[k,j]
        j <- j + 1L
      }
    }

    newbeta <- oldbeta - solve(a = dscore1 + dscore2, 
                               b = colSums(x = score1) + colSums(x = score2))

    error <- sum(abs(x = newbeta-oldbeta))
    iter <- iter + 1L
    oldbeta <- newbeta

  }

  # {nVac x nEventsVac}
  iex <- indVac*EXbeta
 
  inf2 <- matrix(data = 0.0, nrow = nVac, ncol = np+m)
  inf2[DeltaVac == 1L,] <- score2

  k <- 1L
  while (k <= np) {
    # {nVac x nEventsVac}
    tempk <- iex*X[vac,k]

    # {nEventsVac}
    nomk <- colSums(x = tempk)

    inf2[,k] <- inf2[,k] - tempk %*% {1.0/denom} + iex %*% {nomk/{denom^2}}
    k <- k + 1L
  }

  for (k in {np+1L}:{np+m}) {
    # {nVac x nEventsVac}
    tempk <- iex*{loc==k}

    # {nEventsVac}
    nomk <- colSums(x = tempk)
    inf2[,k] <- inf2[,k] - tempk %*% {1.0/denom} + iex %*% {nomk/{denom^2}}
  }

  # {n x {np+m}}
  infbeta <- score1
  infbeta[vac,] <- infbeta[vac,] + inf2
  infbeta <- -infbeta %*% solve(a = dscore1 + dscore2)
 
  # {np+m}
  betase <- sqrt(x = diag(x = crossprod(x = infbeta)))

  denom <- rep(1.0, times = nVac)
  denom[DeltaVac == 1L] <- colSums(x = iex)

  # this is problematic for reals
  testpt <- c(0.0, sort(x = unique(x = diseaseTimeAfterVac[DeltaVac==1L])))

  # {m} each {nEventsVac}
  ttt <- matrix(data = 0.0, nrow = nEventsVac, ncol = m)
  for (k in {np+1L}:{np+m}) ttt[,k-np] <- colSums(iex*{loc==k})

  dd1 <- outer(X = diseaseTimeAfterVac, Y = testpt, FUN = "<=")
  dd2 <- outer(X = eventTimesAfterVac, Y = testpt, FUN = "<=")

  VC <- colSums(DeltaVac*dd1/denom)

  Vt <- VC

  tmp <- dd2/{denom[DeltaVac == 1L]^2}

  temp <- dd1*{DeltaVac/denom} - iex %*% tmp

  SVC <- matrix(data = 0.0, nrow = n, ncol = ncol(x = temp))
  SVC[vac,] <- temp

  if (np > 0L) {
    # {np} each {nEventsVac}
    uuu <- apply(X = X[vac,,drop=FALSE], 
                 MARGIN = 2L,
                 FUN = function(x,y) {colSums(y*x)},
                 y = iex)
    tmp2 <- tcrossprod(x = crossprod(x = tmp, y = uuu), y = infbeta[,1L:np])
  } else {
    tmp2 <- 0.0
  }
  tmp3 <- tcrossprod(x = crossprod(x = tmp, y = ttt), 
                     y = infbeta[,{np+1L}:{np+m}])

  SVC <- SVC - t(x = tmp2 + tmp3)

  sdVC <- sqrt(x = colSums(x = SVC^2))

  ind <- VC > 1e-8
  VC[ind] <- VC[ind] / testpt[ind]
  sdVC[ind] <- sdVC[ind] / testpt[ind]
  logsdVC <- numeric(length = length(x = testpt))
  logsdVC[ind] <- sdVC[ind] / VC[ind]

  outputData <- cbind(testpt, 1.0-VC, sdVC, 
                      1.0-VC*exp(x = -1.96*logsdVC),  
                      1.0-VC*exp(x = 1.96*logsdVC))
  
  colnames(x = outputData) <- c("time", "VE_a", "se", "lower .95", "upper .95")
  
  rr_Combined <- {1.0 - outputData[,2L]}*outputData[,1L]

  outputData <- cbind(outputData, 
                      "hazardRatio" = c(0.0, diff(x = rr_Combined)))

  # remove zero if present
  timePts <- sort(x = unique(x = c(0.0, timePts)))
  timePts <- timePts[-1L]
  nt <- length(x = timePts)

  VE_ave <- numeric(length = nt)
  VE_ave_up <- numeric(length = nt)
  VE_ave_low <- numeric(length = nt)
  diff_se <- numeric(length = nt)

  lowInd <- 1L
  lowTime <- 0.0

  for (k in 1L:nt) {
    # number of time points at or below the upper bound
    ind <- sum(outputData[,1L] <= timePts[k])

    # size of interval
    dt <- timePts[k] - lowTime

    # difference between vaccine efficacy at the two bounds
    diff <- Vt[ind] - Vt[lowInd]

    # standard error
    tmp <- SVC[,ind] - SVC[,lowInd]
    diff_se[k] <- sqrt(x = sum(tmp^2))

    if (diff > 0.0) {
      VE_ave[k] <- 1.0 - diff / dt
      VE_ave_up[k] <- 1.0 - diff*exp(x = -1.96*diff_se[k]/diff)/dt
      VE_ave_low[k] <- 1.0 - diff*exp(x = 1.96*diff_se[k]/diff)/dt
    }

    diff_se[k] <- diff_se[k] / dt

    lowInd <- ind
    lowTime <- timePts[k]
  }

  outputInterval <- cbind("left" = c(0.0,timePts[-nt]), 
                          "right" = timePts,
                          "VE_a" = VE_ave,
                          "se" = diff_se,
                          "lower .95" = VE_ave_low,
                          "upper .95" = VE_ave_up)
    
  rownames(x = outputInterval) <- NULL

  if (np > 0L) {
    
    beta <- newbeta[1L:np]
    
    # output coefficient
    # beta, exp(beta), se, z, p-value, lower 95% CI, upper 95% CI
    xBeta <- beta / betase[1L:np]
    
    outputBeta <- cbind(beta, 
                        betase[1L:np], 
                        xBeta, 
                        2.0*{1.0 - stats::pnorm(q = abs(x = xBeta), 
                                                mean = 0.0, sd = 1.0)},
                        exp(x = beta), 
                        exp(x = beta - 1.96*betase[1L:np]), 
                        exp(x = beta + 1.96*betase[1L:np]))
    
    colnames(x = outputBeta) <- c("coef",
                                  "se(coef)", "z", "Pr(>|z|)",
                                  "exp(coef)",
                                  "lower .95", "upper .95")
    
    rownames(x = outputBeta) <- colnames(x = X)
    
    
  } else {
    outputBeta <- NULL
  }

  return( list("covariates" = outputBeta,
               "vaccine" = list("efficacy" = outputData, 
                                "period_efficacy" = outputInterval)) )
}

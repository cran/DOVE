# note to developer: beta is a matrix
.postVE <- function(X, 
                    beta, 
                    loc, 
                    iex, 
                    DeltaVac, 
                    diseaseTimeAfterVac, 
                    eventTimesAfterVac, 
                    vac) {

  # number of covariates
  np <- ifelse(test = is.null(x = X), yes = 0L, no = ncol(x = X))

  m <- ncol(x = beta) - np
  n <- length(x = vac)

  # number of vaccinated participants
  nVac <- nrow(x = iex)

  # number of events experienced by vaccinated participants
  nEventsVac <- ncol(x = iex)

  denom <- rep(1.0, times = nVac)
  denom[DeltaVac == 1L] <- colSums(x = iex)

  # this could be problematic for reals
  testpt <- sort(x = unique(x = diseaseTimeAfterVac[DeltaVac==1L]))
  if (testpt[1L] > 1e-8) testpt <- c(0.0, testpt)

  # {m} each {nEventsVac}
  ttt <- matrix(data = 0.0, nrow = nEventsVac, ncol = m)
  for (k in {np+1L}:{np+m}) ttt[,k-np] <- colSums(iex*{loc==k})

  dd1 <- outer(X = diseaseTimeAfterVac, Y = testpt+1e-8, FUN = "<")
  dd2 <- outer(X = eventTimesAfterVac, Y = testpt+1e-8, FUN = "<")

  VC <- colSums(x = DeltaVac*dd1/denom)

  Vt <- VC

  tmp <- dd2 / {denom[DeltaVac == 1L]^2}

  temp <- dd1*{DeltaVac/denom} - iex %*% tmp

  SVC <- matrix(data = 0.0, nrow = n, ncol = ncol(x = temp))
  SVC[vac,] <- temp

  if (np > 0L) {
    # {np} each {nEventsVac}
    uuu <- apply(X = X[vac,,drop=FALSE], 
                 MARGIN = 2L,
                 FUN = function(x,y) {colSums(x = y*x)},
                 y = iex)
    tmp2 <- tcrossprod(x = crossprod(x = tmp, y = uuu), y = beta[,1L:np])
  } else {
    tmp2 <- 0.0
  }
  tmp3 <- tcrossprod(x = crossprod(x = tmp, y = ttt), 
                     y = beta[,{np+1L}:{np+m}])

  SVC <- SVC - t(x = tmp2 + tmp3)

  sdVC <- sqrt(x = colSums(x = SVC^2))

  ind <- VC > 1e-8
  VC[ind] <- VC[ind] / testpt[ind]
  sdVC[ind] <- sdVC[ind] / testpt[ind]
  logsdVC <- numeric(length = length(x = testpt))
  logsdVC[ind] <- sdVC[ind] / VC[ind]

  outputData <- cbind(testpt, 1.0-VC, sdVC, 
                      1.0-VC*exp(x = 1.96*logsdVC),  
                      1.0-VC*exp(x = -1.96*logsdVC))
  
  # 7/30/21: restrict VE_a to be 0 at day 0
  outputData[1L,c(2L,4L,5L)] <- 0.0
  
  colnames(x = outputData) <- c("time", "VE_a", "se", "lower .95", "upper .95")
  
  rr_Combined <- {1.0 - outputData[,2L]}*outputData[,1L]

  # 7/30/21: change "hazardRatio" at day 0 from 0.0 to 1.0
  outputData <- cbind(outputData, 
                      "hazardRatio" = c(1.0, diff(x = rr_Combined)))

  return( list("VE" = outputData, "Vt" = Vt, "SVC" = SVC) )

}


.postInterval <- function(timePts, times, VE, SVC) {

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
    ind <- sum(times < {timePts[k]+1e-8})

    # size of interval
    dt <- timePts[k] - lowTime

    # difference between vaccine efficacy at the two bounds
    diff <- VE[ind] - VE[lowInd]

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

  return( outputInterval )

}

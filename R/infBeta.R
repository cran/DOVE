# @param X A matrix object or NULL. The covariates for the vaccinated
#          participants. {nVac x np}
# @param beta A vector object. The estimated parameters. {np+m}
# @param loc A matrix object. Indicator of which component of beta to
#          to use for each bin. {nVac x nEventsVac}
# @param iex A matrix object.  {nVac x nEventsVac}
# @param DeltaVac A vector object. Indicator of a vaccinated participant
#          experienced an event.
# @param vac A vector object. Indicator of a vaccinated participant.
# @param s1 Value object returned by .score1()
# @param s2 Value object returned by .score2()
setGeneric(name = ".infBeta",
           def = function(X, ...) { standardGeneric(".infBeta") })

setMethod(f = ".infBeta",
          signature = c(X = "ANY"),
          definition = function(X, ...) { 
              stop("input type not supported") 
            })

setMethod(f = ".infBeta",
          signature = c(X = "matrix"),
          definition = function(X, loc, iex, DeltaVac, vac, s1, s2) {

              np <- ncol(x = X)
              m <- ncol(x = s2$score) - np

              infCov <- .infBeta_cov(X = X[vac,,drop=FALSE], iex = iex)

              infm <- .infBeta_m(iex = iex, loc = loc-np, m = m)

              inf2 <- cbind(infCov, infm)
              inf2[DeltaVac == 1L,] <- inf2[DeltaVac == 1L,] + s2$score

              # {n x {np+m}}
              infbeta <- s1$score
              infbeta[vac,] <- infbeta[vac,] + inf2

              inv <- tryCatch(expr = solve(a = s1$dscore + s2$dscore),
                              error = function(e) {
                                        message("unable to invert information matrix")
                                        stop(e$message, call. = FALSE)
                                       })

              infbeta <- -infbeta %*% inv

              return( infbeta )
            })

setMethod(f = ".infBeta",
          signature = c(X = "NULL"),
          definition = function(X, loc, iex, DeltaVac, vac, s1, s2) {

              m <- ncol(x = s2$score)

              inf2 <- .infBeta_m(iex = iex, loc = loc, m = m)

              inf2[DeltaVac == 1L,] <- inf2[DeltaVac == 1L,] + s2$score

              # {n x m}
              infbeta <- s1$score
              infbeta[vac,] <- infbeta[vac,] + inf2

              inv <- tryCatch(expr = solve(a = s1$dscore + s2$dscore),
                              error = function(e) {
                                        message("unable to invert information matrix")
                                        stop(e$message, call. = FALSE)
                                       })

              infbeta <- -infbeta %*% inv

              return( infbeta )
            })

.infBeta_cov <- function(X, iex) {

  # number of covariates
  np <- ncol(x = X)

  # number of vaccinated participants
  nVac <- nrow(x = iex)

  # {nEventsVac}
  denom <- 1.0 / colSums(x = iex)

  res <- matrix(data = 0.0, nrow = nVac, ncol = np)

  for (k in 1L:np) {

    # {nVac x nEventsVac}
    tempk <- iex*X[,k]

    # {nEventsVac}
    nomk <- colSums(x = tempk)

    res[,k] <- - tempk %*% denom + iex %*% {nomk * {denom^2}}
  }

  return( res )

}

#loc takes values 0:m
.infBeta_m <- function(iex, loc, m) {

  denom <- 1.0 / colSums(x = iex)
  nVac <- nrow(x = iex)

  res <- matrix(data = 0.0, nrow = nVac, ncol = m)

  for (k in 1L:m) {
 
   # {nVac x nEventsVac}
    tempk <- iex*{loc == k}

    # {nEventsVac}
    nomk <- colSums(x = tempk)

    res[,k] <- - tempk %*% denom + iex %*% {nomk * {denom^2}}

  }

  return( res )
}

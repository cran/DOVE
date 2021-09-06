#' @include infBeta.R postBeta.R postInterval.R
setGeneric(name = ".postProcessDove",
           def = function(X, ...) { standardGeneric(".postProcessDove") })

setMethod(f = ".postProcessDove",
          signature = c(X = "ANY"),
          definition = function(X, ...) { 
              stop("input type not supported") 
            })

setMethod(f = ".postProcessDove",
          signature = c(X = "matrix"),
          definition = function(X, 
                                beta, 
                                loc, 
                                DeltaVac, 
                                s1, 
                                s2, 
                                diseaseTimeAfterVac, 
                                eventTimesAfterVac, 
                                vac,
                                timePts,
                                indVac) {

              # number of covariates
              np <- ncol(x = X)

              # number of vaccinated participants
              nVac <- nrow(x = loc)

              # number of events experienced by vaccinated participants
              nEventsVac <- ncol(x = loc)

              # {n}
              zbeta <- drop(x = X[vac,,drop = FALSE] %*% beta[1L:np])

              # {nVac x nEventsVac}
              tmp <- matrix(data = 0L, nrow = nVac, ncol = nEventsVac)
              tmp[loc > np] <- beta[loc[loc > np]]

              # {nVac x nEventsVac}
              EXbeta <- exp(x = zbeta + tmp)

              # {nVac x nEventsVac}
              iex <- indVac*EXbeta

              return( .step(X = X,
                            beta = beta,
                            iex = iex, 
                            loc = loc, 
                            DeltaVac = DeltaVac, 
                            s1 = s1, 
                            s2 = s2, 
                            diseaseTimeAfterVac = diseaseTimeAfterVac, 
                            eventTimesAfterVac = eventTimesAfterVac, 
                            vac = vac,
                            timePts = timePts) )
            })

setMethod(f = ".postProcessDove",
          signature = c(X = "NULL"),
          definition = function(X, 
                                beta, 
                                loc, 
                                DeltaVac, 
                                s1, 
                                s2, 
                                diseaseTimeAfterVac, 
                                eventTimesAfterVac, 
                                vac,
                                timePts,
                                indVac) {

              # number of vaccinated participants
              nVac <- nrow(x = loc)

              # number of events experienced by vaccinated participants
              nEventsVac <- ncol(x = loc)

              # {nVac x nEventsVac}
              EXbeta <- matrix(data = 1.0, nrow = nVac, ncol = nEventsVac)
              EXbeta[loc > 0L] <- exp(x = beta[loc[loc > 0]])

              # {nVac x nEventsVac}
              iex <- indVac*EXbeta

              return( .step(X = X,
                            beta = beta,
                            iex = iex, 
                            loc = loc, 
                            DeltaVac = DeltaVac, 
                            s1 = s1, 
                            s2 = s2, 
                            diseaseTimeAfterVac = diseaseTimeAfterVac, 
                            eventTimesAfterVac = eventTimesAfterVac, 
                            vac = vac,
                            timePts = timePts) )

            })

#' @include infBeta.R postBeta.R postInterval.R postBeta.R
.step <- function(X,
                  beta,
                  iex, 
                  loc, 
                  DeltaVac, 
                  s1, 
                  s2, 
                  diseaseTimeAfterVac, 
                  eventTimesAfterVac, 
                  vac,
                  timePts) {

  infbeta <- .infBeta(X = X,
                      loc = loc,
                      iex = iex,
                      DeltaVac = DeltaVac,
                      vac = vac,
                      s1 = s1,
                      s2 = s2)

  outputVE <- .postVE(X = X,
                      beta = infbeta,
                      loc = loc,
                      iex = iex,
                      DeltaVac = DeltaVac,
                      diseaseTimeAfterVac = diseaseTimeAfterVac,
                      eventTimesAfterVac = eventTimesAfterVac,
                      vac = vac)

  outputInterval <- .postInterval(timePts = timePts,
                                  times = outputVE$VE[,1L],
                                  VE = outputVE$Vt,
                                  SVC = outputVE$SVC)

  if (!is.null(x = X)) {
    nms <- colnames(x = X)
    if (length(x = nms) == 0L) nms <- paste0("x",1L:ncol(x = X))
  } else {
    nms <- NULL
  }

  outputBeta <- .postBeta(infbeta = infbeta, beta = beta, names = nms)

  return( list("covariates" = outputBeta,
               "vaccine" = list("efficacy" = outputVE$VE, 
                                "period_efficacy" = outputInterval)) )
}

.iteration <- function(score1Prep,
                       X, 
                       indVac, 
                       loc, 
                       score2Prep,
                       m, 
                       vac) {

  np <- ifelse(test = is.null(x = X), yes = 0L, no = ncol(x = X))
  n <- length(x = vac)
  nEventsVac <- ncol(x = indVac)
  nVac <- sum(vac)

  oldbeta <- numeric(length = np+m)
  iter <- 0L 
  error <- 1.0 
  maxiter <- 500L
  epsilon <- 0.0001

  dscore1 <- matrix(data = 0.0, nrow = np+m, ncol = np+m)
  score1 <- matrix(data = 0.0, nrow = n, ncol = np+m)

  iem <- {np+1L}:{np+m}

  while ({error > epsilon} && {iter < maxiter}) {

    # {n}
    if (np > 0L) {
      zbeta <- drop(x = X %*% oldbeta[1L:np])
    } else {
      zbeta <- numeric(length = n)
    }

    # {n x m}
    EXbeta <- exp(x = outer(X = zbeta, Y = oldbeta[iem], FUN = "+"))

    # {n x m}
    exbt <- EXbeta*score1Prep$dt

    score1 <- score1Prep$DZBBy - cbind(rowSums(exbt)*X,exbt) 

    if (np > 0L) {
      for (k in 1L:np) {
        # {n x m}
        tmp <- X[,k]*exbt

        for (j in k:np) {
          dscore1[k,j] <- -sum(X[,j]*tmp)
          dscore1[j,k] <- dscore1[k,j]
        }
        for (j in iem) {
          dscore1[k,j] <- -sum(tmp[,j-np])
          dscore1[j,k] <- dscore1[k,j]
        }
      }
    }

    dscore1[cbind(iem,iem)] <- -colSums(x = exbt)

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

    iexb <- sweep(x = iexb, MARGIN = 2L, STATS = denom, FUN = "/")

    # {nEventsVac x np}
    if (np > 0L) {
      uuu <- crossprod(x = iexb, y = X[vac,,drop = FALSE])
    }

    ttt <- matrix(data = 0.0, nrow = nEventsVac, ncol = m)
    for (i in 1L:m) {
      for (j in 1L:nEventsVac) {
        use <- loc[,j] == {np+i}
        ttt[j,i] = sum(iexb[use,j])
      }
    }

    score2 <- score2Prep - cbind(uuu, ttt)

    if (np > 0L) {

      dscore2[1L:np,] <- crossprod(uuu, cbind(uuu,ttt))
      dscore2[,1L:np] <- t(x = dscore2[1L:np,])

      for (k in 1L:np) {

        # d exp(Z beta + gamma) / d beta_k
        # {nVac x nEventsVac}
        tempk <- iexb*X[vac,k]

        # {nEventsVac x np}
        nomjk <- colSums(x = crossprod(x = tempk, y = X[vac,,drop=FALSE]))
        dscore2[k,k:np] <- dscore2[k,k:np] - nomjk[k:np]
        dscore2[k:np,k] <- dscore2[k,k:np]

        for (j in {np+1L}:{np+m}) {
          # {nEventsVac}
          dnomjk <- sum(tempk[loc == j])

          dscore2[k,j] <- dscore2[k,j] - dnomjk
          dscore2[j,k] <- dscore2[k,j]
        }
      }
    }

    dscore2[cbind(iem,iem)] <- -colSums(x = ttt)
    dscore2[iem,iem] <- dscore2[iem,iem] + crossprod(x = ttt)

    inv <- tryCatch(expr = solve(a = dscore1 + dscore2, 
                                 b = colSums(x = score1) + colSums(x = score2)),
                    error = function(e) {
                              message("unable to invert information matrix")
                              stop(e$message, call. = FALSE)
                             })

    newbeta <- oldbeta - inv

    error <- sum(abs(x = newbeta-oldbeta))
    iter <- iter + 1L
    oldbeta <- newbeta

  }

  message("method converged after ", iter, " iterations")


  return( list(oldbeta, score1, dscore1, score2, dscore2) )
}

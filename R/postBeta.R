.postBeta <- function(infbeta, beta, names) {

  if (!is.null(x = names)) {
    np <- length(x = names)

    betase <- sqrt(x = diag(x = crossprod(x = infbeta)))
   
    # output coefficient
    # beta, exp(beta), se, z, p-value, lower 95% CI, upper 95% CI
    xBeta <- beta[1L:np] / betase[1L:np]
    
    outputBeta <- cbind(beta[1L:np], 
                        betase[1L:np], 
                        xBeta, 
                        2.0*{1.0 - stats::pnorm(q = abs(x = xBeta), 
                                                mean = 0.0, sd = 1.0)},
                        exp(x = beta[1L:np]), 
                        exp(x = beta[1L:np] - 1.96*betase[1L:np]), 
                        exp(x = beta[1L:np] + 1.96*betase[1L:np]))
    
    colnames(x = outputBeta) <- c("coef",
                                  "se(coef)", "z", "Pr(>|z|)",
                                  "exp(coef)",
                                  "lower .95", "upper .95")
    
    rownames(x = outputBeta) <- names
    
    
  } else {
    outputBeta <- NA
  }

  return( outputBeta )

}

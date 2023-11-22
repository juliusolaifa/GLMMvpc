vpc.lrt <- function(vpcObj, nullfixed, nullrandom) {
  logLik1 <- vpcObj$logLikelihood
  m0 <- GLMMadaptive::mixed_model(fixed = nullfixed, random = nullrandom,
                                       data = vpcObj$modObj$data, 
                                       family = vpcObj$modObj$family)
  logLik0 <- m0$logLik
  lambda <- -2*(logLik0 - logLik1)
  pval <- p_value()

  p_value <- function(bou_in, nob_in, nob_nu) {
    if(bou_in == 1 && nob_in == 0 && nob_nu == 2) {
      0.5*stats::pchisq(lambda, df = 1, lower.tail = FALSE) + 
        0.5*(lambda == 0)
    }else if(bou_in == 1 && nob_in == 2 && nob_nu == 3) {
      0.5*stats::pchisq(lambda, df = 1, lower.tail = FALSE) +
        0.5*stats::pchisq(lambda, df = 2, lower.tail = FALSE)
    }else if(bou_in == 1 && nob_in == 3 && nob_nu == 2) {
      0.5*stats::pchisq(lambda, df = 1, lower.tail = FALSE) +
        0.5*stats::pchisq(lambda, df = 2, lower.tail = FALSE)
    }else if(bou_in == 2 && nob_in == 1 && nob_nu == 3) {
      V <- stats::vcov(vpcObj$modObj)
      vbvar <- matrix(c(V[3,3],V[3,5],
                        V[5,3],V[5,5]),
                      nrow = 2, ncol = 2, byrow = TRUE)
      pi1 <- as.numeric(mvtnorm::pmvnorm(lower = c(0,0), upper = c(Inf,Inf),
                                         mean = c(0,0), sigma = vbvar))
      pi4 <- 0.5 - pi1
      pi1*stats::pchisq(lambda, df = 1, lower.tail = FALSE) +
        0.5*stats::pchisq(lambda, df = 2, lower.tail = FALSE) +
        pi4*stats::pchisq(lambda, df = 3, lower.tail = FALSE) 
    }else{
      stop("`vpc.lrt` not implemented for this case")
    }
  }
}
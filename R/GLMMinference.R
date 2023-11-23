lrtest <- function(x, ...) {
  UseMethod("lrtest")
}


lrtest.vpcscore <- function(vpcObj, H0fixed, H0random = NULL, alpha = 0.05) {
  logLik1 <- vpcObj$logLikelihood
  
  if (is.null(H0random)) {
    family <- match_family(vpcObj$modObj$family$family)
    m0 <- suppressWarnings(glmmTMB::glmmTMB(formula = H0fixed, family = family, 
              data = vpcObj$modObj$data))
    m0$nvarcomp_bounded <- 0
    m0$nfixed <- length(m0$obj$par)
 
  } else{
    m0 <- glmm.fit(fixed = H0fixed, random = H0random,
                                    data = vpcObj$modObj$data, 
                                    family = vpcObj$modObj$family)
  }
  
  logLik0 <- as.numeric(stats::logLik(m0))
  lambda <- -2*(logLik0 - logLik1)

    if(is.null(H0random) && vpcObj$nvarcomp_bounded == 1) {
      pval <- 0.5*stats::pchisq(lambda, df = 1, lower.tail = FALSE) + 
        0.5*(lambda == 0)
    }else if(vpcObj$nvarcomp_bounded == 2 && m0$nvarcomp_bounded == 1) {
      if(vpcObj$nfixef == m0$nfixef){
        pval <- 0.5*stats::pchisq(lambda, df = 1, lower.tail = FALSE) +
        0.5*stats::pchisq(lambda, df = 2, lower.tail = FALSE)
      }else if(vpcObj$nfixef == m0$nfixef + 1) {
        pval <- 0.5*stats::pchisq(lambda, df = 2, lower.tail = FALSE) +
        0.5*stats::pchisq(lambda, df = 3, lower.tail = FALSE)
      }
    }else if(is.null(H0random) && vpcObj$nvarcomp_bounded == 2) {
      V <- stats::vcov(vpcObj$modObj)
      vbvar <- matrix(c(V[3,3],V[3,5],
                        V[5,3],V[5,5]),
                      nrow = 2, ncol = 2, byrow = TRUE)
      pi1 <- as.numeric(mvtnorm::pmvnorm(lower = c(0,0), upper = c(Inf,Inf),
                                         mean = c(0,0), sigma = vbvar))
      pi4 <- 0.5 - pi1
      pval <- pi1*stats::pchisq(lambda, df = 1, lower.tail = FALSE) +
        0.5*stats::pchisq(lambda, df = 2, lower.tail = FALSE) +
        pi4*stats::pchisq(lambda, df = 3, lower.tail = FALSE) 
    }
    result <- list("lambda" = lambda, "p_value" = pval, "alpha" = alpha)
    class(result) <- "glmmlrtest"
    return(result)
}

print.glmmlrtest <- function(testObj) {
  cat("Likelihood Ratio Test\n")
  cat("---------------------\n")
  cat("Lambda: ")
  cat(testObj$lambda)
  cat("\n")
  cat("P-value: ")
  cat(testObj$p_value)
  cat("\n")
  if(testObj$p_value < testObj$alpha)
    cat("Reject the Null Hypothesis")
  else
    cat("Do not reject the Null hypothesis")
}

match_family <- function(fam) {
  adaptive_fam <- c("negative binomial")
  glmmTMB_fam <- c("nbinom2")
  index <- match(fam, adaptive_fam)
  if (!is.na(index)) {
    return(glmmTMB_fam[index])
  } else {
    return(NA)
  }
}


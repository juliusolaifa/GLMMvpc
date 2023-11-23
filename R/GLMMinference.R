lrtest <- function(vpcObj, ...) {
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

############################################
confint.vpcScore <- function(vpcObj, level=0.95) {
  
  wrapper_function <- function(params, fitmod, X) {
    len_beta <- length(fitmod$beta)
    len_Sigma <- choose(nrow(fitmod$Sigma), 2) + nrow(fitmod$Sigma)
    
    beta <- params[1:len_beta]
    sigma_vector <- params[(len_beta + 1):(len_beta + len_Sigma)]
    phi <- params[(len_beta + len_Sigma + 1):length(params)]
    
    Sigma <- matrix(0, nrow = nrow(fitmod$Sigma), ncol = nrow(fitmod$Sigma))
    Sigma[upper.tri(Sigma, diag = TRUE)] <- sigma_vector
    Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
    
    family <- fitmod$family
    link <- fitmod$link
    p <- fitmod$p
    # Adjust this call according to the actual parameters expected by glmm.vpc
    return(glmm.vpc(X, beta = beta, Sigma = Sigma, phi = phi, family = family, link = link, p = p))
  }
  fitmod <- vpcObj$fitmod
  sigma_vector <- fitmod$Sigma[upper.tri(fitmod$Sigma, diag = TRUE)]
  param_vector <- c(fitmod$beta, sigma_vector, fitmod$phi)
  J <- numDeriv::grad(func = function(p) wrapper_function(p, fitmod, vpcObj$X), x = param_vector)
  V_n <- vcov(fitmod) / fitmod$n 
  se <- J^T%*%V_n%*%J
  alpha_2 <- (1 - level) / 2 
  z_alpha <- -qnorm(alpha_2)
  lower <- vpcObj$vpc - z_alpha*se
  upper <- vpcObj$vpc + z_alpha*se
  
  result <- list("lower" = lower, "pointestimate" = vpcObj$vpc, "upper" = upper, "level" = level)
  class(result) <- "vpcConfint"
  result
}

print.vpcConfint <- function(vpcObj) {
  cat("Confidence Interval\n")
  cat("-------------------\n")
  cat("lower: ")
  cat(vpcObj$lower)
  cat("\n")
  cat("point estimate: ")
  cat(vpcObj$pointestimate)
  cat("\n")
  cat("upper: ")
  cat(vpcObj$upper)
  cat("\n")
  cat("level: ")
  cat(vpcObj$level)
}


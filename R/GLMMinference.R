lrtest <- function(vpcObj, ...) {
  UseMethod("lrtest")
}


lrtest.vpcScore <- function(vpcObj, H0fixed, H0random = NULL, alpha = 0.05) {
  
  if (is.null(H0random)) {
    family <- match_family(vpcObj$fitmod$family$family)
    m0 <- suppressWarnings(glmmTMB::glmmTMB(formula = H0fixed, family = family, 
              data = vpcObj$fitmod$modObj$data))
    m0$nvarcomp_bounded <- 0
    m0$nfixed <- length(m0$obj$par)
 
  } else{
    family <- gsub(" ", ".", vpcObj$fitmod$family$family)
    m0 <- glmm.fit(fixed = H0fixed, random = H0random,
                                    data = vpcObj$fitmod$modObj$data, 
                                    family = family)
  }
  
  logLik1 <- as.numeric(stats::logLik(vpcObj$fitmod))
  logLik0 <- as.numeric(stats::logLik(m0))
  lambda <- -2*(logLik0 - logLik1)

    if(is.null(H0random) && vpcObj$nvarcomp_bounded == 1) {
      pval <- 0.5*stats::pchisq(lambda, df = 1, lower.tail = FALSE) + 
        0.5*(lambda == 0)
    }else if(vpcObj$fitmod$nvarcomp_bounded == 2 && m0$nvarcomp_bounded == 1) {
      if(vpcObj$fitmod$nfixef == m0$nfixef){
        pval <- 0.5*stats::pchisq(lambda, df = 1, lower.tail = FALSE) +
        0.5*stats::pchisq(lambda, df = 2, lower.tail = FALSE)
      }else if(vpcObj$fitmod$nfixef == m0$nfixef + 1) {
        pval <- 0.5*stats::pchisq(lambda, df = 2, lower.tail = FALSE) +
        0.5*stats::pchisq(lambda, df = 3, lower.tail = FALSE)
      }
    }else if(is.null(H0random) && vpcObj$fitmod$nvarcomp_bounded == 2) {
      V <- stats::vcov(vpcObj$fitmod)
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
    print("Reject the Null Hypothesis")
  else
    print("Do not reject the Null hypothesis")
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
  
  #This function is needed to produce the jacobian 
  # of vpc at the values of the fitted parameters
  wrapper_function <- function(params, fitmod, X) {
    
    nbeta <- length(fitmod$beta)
    #obtain number of unique variance component
    nSigma <- sum(choose(fitmod$nvarcomp_bounded, 2) , fitmod$nvarcomp_bounded)
    beta <- fitmod$beta
    sigma_vector <- params[sum(nbeta, 1):sum(nbeta, nSigma)]
    phi <- fitmod$phi
    
    Sigma <- matrix(0, nrow = nrow(fitmod$Sigma), ncol = nrow(fitmod$Sigma))
    Sigma[upper.tri(Sigma, diag = TRUE)] <- sigma_vector
    Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
    
    family <- fitmod$family$family
    link <- fitmod$link
    p <- fitmod$p
    return(glmm.vpc(X, beta = beta, Sigma = Sigma, phi = phi, family = family, link = link, p = p))
  }
  
  fitmod <- vpcObj$fitmod
  sigma_vector <- fitmod$Sigma[upper.tri(fitmod$Sigma, diag = TRUE)]
  param_vector <- c(fitmod$beta, sigma_vector, fitmod$phi)
  
  J <- numDeriv::grad(func = function(par) wrapper_function(par, fitmod, vpcObj$X), x = param_vector)
  V_n <- f.theta(fitmod)$Sig / fitmod$n 
  se <- J^T%*%V_n%*%J
  alpha_2 <- (1 - level) / 2 
  z_alpha <- -stats::qnorm(alpha_2)
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


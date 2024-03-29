glmm <- function(formula, data, family) {
    modObj <- glmmTMB::glmmTMB(formula=formula, data=data, family=family)
    beta <- unname(fixef(modObj)$cond)
    Sigma <- unname(unlist(VarCorr(modObj)$cond))
    phi <- glmmTMB::sigma(modObj)
    p <- unname(glmmTMB::family_params(modObj))
    family <- modObj$modelInfo$family
    vcov  <-  unname(stats::vcov(modObj, full = TRUE))
    logLik <- as.numeric(logLik(modObj))
    group_var <- modObj$modelInfo$grpVar
    nfixed <- length(c(beta,phi))
    nvarcomp_bounded <- length(attr(VarCorr(modObj)$cond[[group_var]], "stddev"))
    
    result <- list("beta"=beta, "Sigma"=Sigma, "phi"=phi, "p"=p,
                   "family"=family$family, "link"=family$link, 
                   "invlink"=family$linkinv,"nfixed" = nfixed, 
                   "vcov" = vcov,"nvarcomp_bounded"=nvarcomp_bounded, 
                   "n" = nrow(data),  "logLikelihood"=logLik, 
                   "modObj"=modObj)
    class(result) <- "glmmMod"
    result
}

vpc <- function(fitObj, x) {
  # Extract parameters from fitObj
  b0 <- fitObj$beta[1]
  b1 <- if (fitObj$nfixed == 3) fitObj$beta[2] else 0
  phi <- fitObj$phi
  p <- fitObj$p
  n <- fitObj$n
  sig11 <- fitObj$Sigma[1]
  sig12 <- ifelse(fitObj$nvarcomp_bounded == 2, fitObj$Sigma[2], 0)
  sig22 <- ifelse(fitObj$nvarcomp_bounded == 2, fitObj$Sigma[4], 0)
  
  # Define the vpc function
  vpc_func = get(paste0("vpc_",fitObj$family))
  # Compute the value of vpc
  vpc_value <- switch(fitObj$family,
                      "nbinom2" = vpc_func(b0,b1,phi,sig11,sig12,sig22,x),
                      "tweedie" = vpc_func(b0,b1,phi,sig11,sig12,sig22,p,x)
  )
  
  # Wrapper function for gradient calculation
  vpc_wrapper <- function(params) {
    switch(fitObj$family,
           "nbinom2" = vpc_func(params[1], params[2], params[3], params[4], 
                                params[5], params[6], x),
           "tweedie" = vpc_func(params[1], params[2], params[3], params[4], 
                                params[5], params[6], params[7], x)
    )
  }
  
  # Parameters vector for gradient calculation
  params <-switch(fitObj$family,
                  "nbinom2" = c(b0, b1, phi, sig11, sig12, sig22),
                  "tweedie" = c(b0, b1, phi, sig11, sig12, sig22, p),
  )
  
  # Calculate gradient
  gradient <- numDeriv::grad(func = vpc_wrapper, x = params)
  
  # Standard Error
  stderr <- sqrt((gradient %*% fitObj$vcov  %*% gradient)/n)
  
  # Return both vpc value and gradient
  result <- list("vpc" = vpc_value, "gradient" = gradient, "stderr" = stderr, "modObj" = fitObj$modObj)
  class(result) <- "vpcObj"
  return(result)
}

vpc_nbinom2 <- function(b0, b1, phi, sig11, sig12, sig22, x) {
  mu <- b0 + b1 * x
  Sigma <- matrix(c(sig11, sig12, sig12, sig22), 2, 2)
  Z <- c(1, x)
  sig <- Z %*% Sigma %*% Z
  cmp1 <- (exp(sig) - 1) * exp(2 * mu + sig) # varLogNormal
  cmp2 <- exp(mu + sig / 2)  # meanLogNormal
  cmp3 <- (1 / phi) * exp(2 * mu + 2 * sig)  # 1/phi X 2ndMomentLogNormal
  return(cmp1 / (cmp1 + cmp2 + cmp3))  # varfunc = mu +mu^2/phi
}

vpc_tweedie <-  function(b0, b1, phi, sig11, sig12, sig22, p, x) {
  mu <- b0 + b1 * x
  Sigma <- matrix(c(sig11, sig12, sig12, sig22), 2, 2)
  Z <- c(1, x)
  sig <- Z %*% Sigma %*% Z
  cmp1 <- (exp(sig) - 1) * exp(2 * mu + sig) # varLogNormal
  cmp2 <- phi* exp(p*mu + 0.5*p^2*sig) #phi X pthMomentLogNormal
  return(cmp1/(cmp1 + cmp2))  # varfunc = phi*mu^p
}

## Define a function to compute critical value?
confint.vpcObj <- function(vpcObj, level=0.95) {
  z_alph_2 <- -stats::qnorm((1-level)/2)
  lcl <- vpcObj$vpc - z_alph_2 * vpcObj$stderr
  ucl <- vpcObj$vpc + z_alph_2 * vpcObj$stderr
  result <- c(lcl, vpcObj$vpc, ucl)
  names(result) <- c("lower", "point estimate", "upper")
  return(result)
}

## TO BE COMPLETED
lrt.vpcObj <- function(vpcObj, H0, alpha=0.05) {
    full_model <- vpcObj$modObj
    reduced_model <- glmm(formula = H0, data = full_model$frame, 
                                        family=full_model$family)
    chi_sq <- anova(full_model, reduced_model)
    chi_sq_stat <- chi_sq$Chisq[2]
}

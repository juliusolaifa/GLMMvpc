vpc.score <- function(fixed, random, data, family, X, Z=NULL) {
  fitmod <- glmm.fit(fixed=fixed, random=random, data=data, family=family)
  vpc <- glmm.vpc(fitmod, X)
  dimnames(vpc) <- list("", "vpc.score")
  fitmod$score <- vpc
  attr(fitmod, "class") <- "vpcscore"
  fitmod
}

glmm.fit <- function(fixed, random, data, family) {
  
  modObj <- GLMMadaptive::mixed_model(fixed, random, data, family)
  attr(m1$D, "L") <- NULL

  result <- list("beta"=modObj$coefficients, "Sigma"=modObj$D, "phi"=modObj$phi, 
                 "family"=modObj$family$family, "link"=modObj$family$link, 
                 "nfixef"=length(c(modObj$coefficients, modObj$phis)), 
                 "nvarcomp_bounded"=dim(modObj$D)[1], "v_chol" = stats::vcov(modObj),
                 "logLikelihood"=modObj$logLik, "modObj"=modObj)
  
  class(result) <- "glmmMod"
  result
}

glmm.vpc <- function(X, Z=NULL, beta, Sigma, phi, family, link) {
  UseMethod("glmm.vpc")
}

glmm.vpc.default <- function(X,Z=NULL,beta,Sigma,phi,family,link) {
  if(is.null(Z)) {
    Z = X
  }
  if (length(X) != length(beta) || length(Z) != ncol(Sigma)) {
    stop("Dimensions of `X and beta` or `Z and Sigma` are 
         not compatible for multiplication.")
  }
  mu <- X%*%beta
  sigm <- Z%*%Sigma%*%Z
  vpc_compute(mu, sigm, phi, family, link)
}

glmm.vpc.glmmMod <- function(modelObj, X, Z=NULL) {
  beta <- modelObj$beta
  Sigma <- modelObj$Sigma
  phi <- modelObj$phi
  family <- modelObj$family
  link <- modelObj$link
  
  glmm.vpc.default(X, Z, beta, Sigma, phi, family, link)
}

vpc_compute <- function(mu, sigm, phi, family, link, p=NULL) {
  
  if(link == "log") {
    inv.mu <- exp(mu + sigm/2)
    inv.var <- (exp(sigm) - 1)*exp(2*mu + sigm)
  }
  switch(family,
         "negative binomial"={
           inv.mu.p <- exp(2*mu + 4*sigm/2)
           return(inv.var / (inv.var + inv.mu + inv.mu.p/phi))
         },
         tweedie={
           if(is.null(p) || p <= 1 || p >= 2) {
             stop("p must be an integer between 1 and 2")
           }
           inv.mu.p <- exp(p*mu + p^2*sigm/2)
           return(inv.var / (inv.var + phi*inv.mu.p))
         },
         compois={
           return(inv.var / (inv.var + phi*inv.mu))
         },
         genpois={
           return(inv.var / (inv.var + phi*inv.mu)) #phi = exp(eta)
         },
         {
           return(paste("VPC not implemented for family: ", family))
         }
  )
}


vcov.glmmMod <-function(modObj) {
  modObj <- modObj$modObj
  beta <- unname(modObj$coefficients)
  sigma_chol <- chol_transf(modObj$D)
  phis <- unname(modObj$phis)
  theta_chol <-c(beta, sigma_chol, phis)
  
  sigma_f_ind <-length(beta)+1
  sigma_l_ind <-length(beta)+ length(sigma_chol)
  
  chol_tocov <-function(theta) {
    sigma <- theta[sigma_f_ind:sigma_l_ind]
    sigma <- chol_transf(sigma)
    sigma <- sigma[upper.tri(sigma, diag = TRUE)]
    theta[sigma_f_ind:sigma_l_ind] <- sigma
    theta
    }
  J <- numDeriv::jacobian(chol_tocov, theta_chol)
  V_chol <- stats::vcov(modObj)
  V <- J%*%V_chol%*% t(J)
  dimnames(V) <-dimnames(V_chol)
  V
}

## Taken from the GLMMadaptive package
chol_transf <- function (x) {
  if (any(is.na(x) | !is.finite(x)))
    stop("NA or infinite values in 'x'.\n")
  if (is.matrix(x)) {
    k <- nrow(x)
    U <- chol(x)
    U[cbind(1:k, 1:k)] <- log(U[cbind(1:k, 1:k)])
    U[upper.tri(U, TRUE)]
  } else {
    nx <- length(x)
    k <- round((-1 + sqrt(1 + 8 * nx))/2)
    mat <- matrix(0, k, k)
    mat[upper.tri(mat, TRUE)] <- x
    mat[cbind(1:k, 1:k)] <- exp(mat[cbind(1:k, 1:k)])
    res <- crossprod(mat)
    attr(res, "L") <- t(mat)[lower.tri(mat, TRUE)]
    res
  }
}

glmm.fit <- function(fixed, random, data, family) {
  
  modObj <- GLMMadaptive::mixed_model(fixed, random, data, family)
  result <- list("beta"=modObj$coefficients, "Sigma"=modObj$D, "phi"=modObj$phi, 
                 "family"=family, "link"=modObj$family$link, 
                 "nfixef"=length(c(modObj$coefficients, modObj$phis)), 
                 "nvarcomp_bounded"=dim(modObj$D)[1], "v_chol" = stats::vcov(modObj),
                 "logLikelihood"=modObj$logLik, "modObj"=modObj, "n" = nrow(data))
  class(result) <- "glmmMod"
  result
}

vpc.score <- function(fixed, random, data, family, X, Z=NULL) {
  fitmod <- glmm.fit(fixed=fixed, random=random, data=data, family=family)
  vpc <- glmm.vpc(fitmod, X)
  dimnames(vpc) <- list("", "vpc.score")
  result <- list("fitmod" = fitmod, "f.theta" = f.theta(fitmod), "vpc" = vpc,
                 "X" = X)
  class(result) <- "vpcScore"
  result
}

glmm.vpc <- function(X, Z=NULL, beta, Sigma, phi, family, link, p) {
  UseMethod("glmm.vpc")
}

glmm.vpc.default <- function(X,Z=NULL,beta,Sigma,phi,family,link, p) {
  if(is.null(Z)) {
    Z = X
  }
  if (length(X) != length(beta) || length(Z) != ncol(Sigma)) {
    stop("Dimensions of `X and beta` or `Z and Sigma` are 
         not compatible for multiplication.")
  }
  mu <- X%*%beta
  sigm <- Z%*%Sigma%*%Z
  vpc_compute(mu, sigm, phi, family, link, p)
}

glmm.vpc.glmmMod <- function(modelObj, X, Z=NULL) {
  beta <- modelObj$beta
  Sigma <- modelObj$Sigma
  phi <- modelObj$phi
  family <- modelObj$family
  link <- modelObj$link
  p <- modelObj$p
  glmm.vpc.default(X, Z, beta, Sigma, phi, family, link, p)
}

vpc_compute <- function(mu, sigm, phi, family, link, p=NULL) {
  
  if(link == "log") {
    inv.mu <- exp(mu + sigm/2)
    inv.var <- (exp(sigm) - 1)*exp(2*mu + sigm)
  }
  switch(family,
         "negative.binomial"={
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
           stop("VPC not implemented for family: ", family)
         }
  )
}

print.vpcScore <- function(vpcSoreObj) {
  cat("Vpc: ")
  cat(vpcSoreObj$vpc)
  cat("\n")
  cat("Distribution\n")
  print(vpcSoreObj$f.theta)
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

f.theta <- function(modObj) {
  V <- stats::vcov(modObj)
  sigma_f_ind <- length(modObj$beta) + 1
  nvar <- modObj$nvarcomp_bounded
  npar <- dim(V)[1]
  
  vbvar <- matrix(c(V[3,3],V[3,5],
                    V[5,3],V[5,5]),
                   nrow = 2, ncol = 2, byrow = TRUE)
  
  if (nvar == 1) {
    V1 <- V2 <- V
    V2 <- V2[-sigma_f_ind, -sigma_f_ind]
    lower1 <- rep(-Inf,npar)
    lower1[sigma_f_ind] <- 0
    f1 <- tmvtnorm::mtmvnorm(mean = rep(0,npar), sigma = V1, 
                   lower = lower1, upper = rep(Inf,npar))
    mu1 <- f1$tmean
    Sigma1 <- f1$tvar
    pi1 <- 0.5
    f2 <- tmvtnorm::mtmvnorm(mean = rep(0,npar-1), sigma = V2, 
                   lower = rep(-Inf,npar-1), upper = rep(Inf,npar-1))
    mu2 <- f2$tmean
    Sigma2 <- f2$tvar
    pi2 <- 0.5
    mu <- pi1*mu1+pi2*mu2
    Sig <- (pi1*(Sigma1 + mu1 %*%t(mu1))+pi2*(Sigma2 + mu2 %*%t(mu2)))- mu %*%t(mu)
    return(list("mu" = mu, "Sig" = Sig))
    # return(list("pis" = list("pi_1" = 0.5, "pi_2" = 0.5),
    #             "means" = list("mu_1" = mu1, "mu_2" = mu2), 
    #             "Sigma" = list("sig_1" = Sigma1, "sig_2" = Sigma2)))
  } else if(nvar == 2) {
    V1 <- V2 <- V3 <- V4 <- V
    V2 <- V2[-sigma_f_ind, -sigma_f_ind]
    V3 <-V3[-(sigma_f_ind+nvar), -(sigma_f_ind+nvar)]
    V4 <- V4[-c(sigma_f_ind,sigma_f_ind+nvar), -c(sigma_f_ind,sigma_f_ind+nvar)]
    lower1 <- rep(-Inf,npar)
    lower1[c(sigma_f_ind,sigma_f_ind+nvar)] <- 0
    f1 <- tmvtnorm::mtmvnorm(mean = rep(0,npar), sigma = V1, 
                   lower = lower1, upper = rep(Inf,npar))
    pi1 <- as.numeric(mvtnorm::pmvnorm(lower = c(0,0), upper = c(Inf,Inf),
                              mean = c(0,0), sigma = vbvar))
    mu1 <- f1$tmean
    Sigma1 <- f1$tvar
    f2 <- tmvtnorm::mtmvnorm(mean = rep(0,npar-1), sigma = V2, 
                  lower = c(-Inf,-Inf,-Inf,0,-Inf), upper = rep(Inf,npar-1))
    pi2 <- 0.25
    mu2 <- f2$tmean; mu2 <- c(mu2[1:(sigma_f_ind-1)], 0, mu2[sigma_f_ind:length(mu2)])
    Sigma2 <- f2$tvar
    Sigma2 <- rbind(Sigma2[1:(sigma_f_ind-1), ],0, Sigma2[sigma_f_ind:nrow(Sigma2),])
    Sigma2 <- cbind(Sigma2[,1:(sigma_f_ind-1)],0, Sigma2[,sigma_f_ind:ncol(Sigma2)])
    f3 <- tmvtnorm::mtmvnorm(mean = rep(0,npar-1), sigma = V3, 
                   lower = c(-Inf,-Inf,0,-Inf,-Inf), upper = rep(Inf,npar-1))
    pi3 <- 0.25
    mu3 <- f3$tmean; mu3 <- c(mu3[1:(sigma_f_ind+1)], 0, mu3[length(mu3)])
    Sigma3 <- f3$tvar
    Sigma3 <- rbind(Sigma3[1:(sigma_f_ind+1), ],0, Sigma3[nrow(Sigma3),])
    Sigma3 <- cbind(Sigma3[,1:(sigma_f_ind+1)],0, Sigma3[,ncol(Sigma3)])
    f4 <- tmvtnorm::mtmvnorm(mean = rep(0,npar-2), sigma = V4, 
                  lower = rep(-Inf,npar-2), upper = rep(Inf,npar-2))
    pi4 <-0.5 - pi1
    mu4 <- f4$tmean; mu4 <- c(mu4[1:(sigma_f_ind-1)], 0, mu4[sigma_f_ind], 0, 
                              mu4[length(mu4)])
    Sigma4 <- f4$tvar
    Sigma4 <- rbind(Sigma4[1:(sigma_f_ind-1), ],0, Sigma4[sigma_f_ind,], 0, 
                    Sigma4[nrow(Sigma4),])
    Sigma4 <- cbind(Sigma4[,1:(sigma_f_ind-1)],0, Sigma4[,sigma_f_ind], 0, 
                    Sigma4[,ncol(Sigma4)])
    mu <- pi1*mu1+pi2*mu2+pi3*mu3+pi4*mu4
    Sig <- (pi1*(Sigma1 + mu1 %*%t(mu1))+pi2*(Sigma2 + mu2 %*%t(mu2))+
          pi3*(Sigma3 + mu3 %*%t(mu3))+pi4*(Sigma4 + mu4 %*%t(mu4)))- mu %*%t(mu)
    return(list("mu" = mu, "Sig" = Sig))
    # return(list("pis"=list("pi_1" = pi1, "pi_2" = pi2, "pi_3" = pi3, "pi_4" = pi4), 
    #             "means" = list("mu_1" = mu1, "mu_2" = mu2, "mu_3" = mu3, "mu_4" = mu4), 
    #             "Sigma" = list("sig_1" = Sigma1, "sig_2" = Sigma2, "sig_3" = Sigma3, 
    #                            "sig_4" = Sigma4)))
  }
}

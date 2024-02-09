  simulate <- function(iter.size, x, beta, sigma.u, group.sizes, fitfam, 
                                                formula, family,...) {
  data <- glmmdata(iter.size=iter.size, x=x, beta=beta, sigma.u=sigma.u, 
                   group.sizes=group.sizes, family=family,...)
  args <- list(...)
  vpc_true0 <- ifelse(family=="nbinom2", 
                    vpc_nbinom2(beta[1],beta[2],args$theta,sigma.u[1,1],
                                sigma.u[1,2], sigma.u[2,2], 0),
                    vpc_tweedie(beta[1],beta[2],args$theta,sigma.u[1,1],
                                sigma.u[1,2], sigma.u[2,2], args$power, 0))
  vpc_true1 <- ifelse(family=="nbinom2", 
                      vpc_nbinom2(beta[1],beta[2],args$theta,sigma.u[1,1],
                                  sigma.u[1,2], sigma.u[2,2], 1),
                      vpc_tweedie(beta[1],beta[2],args$theta,sigma.u[1,1],
                                  sigma.u[1,2], sigma.u[2,2], args$power, 1))
  fits <- lapply(1:nrow(data), function(i) glmm(formula=formula, data[i], family=fitfam))
  vpc0 <- lapply(fits, function(fit) vpc(fitObj=fit, x=0)) 
  vpc1 <- lapply(fits, function(fit) vpc(fitObj=fit, x=1))
  
  return(list("true0"=vpc_true0, "true1"=vpc_true1, "est0"=vpc0, "est1"=vpc1, "fits"=fits))
}


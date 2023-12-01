set.seed(+234)
# dataGen <- function(iter, ngroups, matnrow, beta, family, ...) {
#   group.sizes <- rbinom(ngroups,10,0.8)
#   x <- rgen01(group.sizes)
#   Sigma <- rpdmat(matnrow)
#   print(theta)
#   data <- glmmdata(iter, x, beta, Sigma, group.sizes, family, ...)
# }


# for generating random positive definite matrix of dimension n
rpdmat <- function(n, a=0, b=1) {
  A <- matrix(stats::runif(n*n, min=a, max=b), nrow = n)
  A <- A%*%t(A)
  return(A + diag(stats::runif(1), n))
}

# for generating binary covariate equally randomly 
# distributed within groups (random effect)
rgen01 <- function(ns) {
  unlist(
    lapply(ns, function(x) 
      sample(rep_len(c(0,1),length.out=x)
  )))
}

group.sizes <- rbinom(10,10,0.8)
beta <- c(7,9)
sigma.u <- rpdmat(2)
x <- rgen01(group.sizes)
theta = 2
phi = 1.8
iter =10
fixed = y ~ x
random = ~ 1 + x |group
red_random = ~ 1 |group
data1 <- glmmdata(iter, x, beta, sigma.u, group.sizes, "nbinom2", theta=theta)
vpc0 <- vpc.score(fixed = fixed, random = random, data=data1[1], 
          family = GLMMadaptive::negative.binomial(), X = c(1,0))
ci0 <- confint(vpc0)
vpc1 <- vpc.score(fixed = fixed, random = random, data=data1[1], 
                  family = GLMMadaptive::negative.binomial(), X = c(1,1))
ci1 <- confint(vpc1)
lrtest(vpc0, H0fixed=fixed, H0random = red_random, alpha = 0.05)

#GLMMadaptive::mixed_model(fixed = fixed, random = random, data=data1[1], 
#                       family = "negative.binomial")

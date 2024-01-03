library(parallel)
library(foreach)
#library(pbapply)
library(ggplot2)
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


group.sizes <- rbinom(10,6,0.9)
beta <- c(2,3)
sigma.u <- rpdmat(2)#matrix(c(3,1,1,2),2)#rpdmat(2)
x <- rgen01(group.sizes)
# theta = 2
phi = 1.8
iter =100
fixed = y ~ x
random = ~ 1 + x |group
red_random = ~ 1 |group
true0 = as.numeric(glmm.vpc(X = c(1,0),beta=beta,Sigma=sigma.u,phi=phi,
               family="negative binomial",link="log", p=NULL))
true0
true1 = as.numeric(glmm.vpc(X = c(1,1),beta=beta,Sigma=sigma.u,phi=phi,
                            family="negative binomial",link="log", p=NULL))
true1
data1 <- glmmdata(iter, x, beta, sigma.u, group.sizes, "nbinom2", theta=theta)
vpc0 <- vpc.score(fixed = fixed, random = random, data=data1[1],
          family = GLMMadaptive::negative.binomial(), X = c(1,0))
ci0 <- confint(vpc0)
vpc1 <- vpc.score(fixed = fixed, random = random, data=data1[1],
                  family = GLMMadaptive::negative.binomial(), X = c(1,1))
ci1 <- confint(vpc1)
lrt0 <- lrtest(vpc0, H0fixed=fixed, H0random = red_random, alpha = 0.05)

# data1
# data1[1]
# vpc0
# vpc1
ci0
ci1
lrt0

### Parallel
cl <- makeCluster(detectCores()-1)
clusterExport(cl, c("data1", "vpc.score"))

# compute_vpcp <- function(i) {
#   vpc.score(fixed = fixed, random = random, data = data1[i], 
#             family = GLMMadaptive::negative.binomial(), X = c(1, 0))
# }
# results <- pbsapply(1:length(data1), compute_vpcp, cl = cl)
results0 <- foreach(i = 1:nrow(data1), .inorder = FALSE, .combine = "rbind", .packages = c("parallel")) %dopar% {
  vpc.score(fixed = fixed, random = random, data=data1[i],
            family = GLMMadaptive::negative.binomial(), X = c(1,0))$vpc
}

results1 <- foreach(i = 1:nrow(data1), .inorder = FALSE, .combine = "rbind", .packages = c("parallel")) %dopar% {
  vpc.score(fixed = fixed, random = random, data=data1[i],
            family = GLMMadaptive::negative.binomial(), X = c(1,1))$vpc
}
stopCluster(cl)

bias0 <- true0 - results0
bias0 <- data.frame(bias0)
colnames(bias0) <- "Bias"

bias1 <- true1 - results1
bias1 <- data.frame(bias1)
colnames(bias1) <- "Bias"

ggplot(bias0, aes(x = Bias)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Bias of Estimator (Density Plot)", x = "Bias")

ggplot(bias1, aes(x = Bias)) +
  geom_density(fill = "red", alpha = 0.5) +
  labs(title = "Bias of Estimator (Density Plot)", x = "Bias")

par(mfrow=c(2,2))
boxplot(bias0, main = "Bias of Estimator (Boxplot)0", ylab = "Bias")
boxplot(bias1, main = "Bias of Estimator (Boxplot)1", ylab = "Bias")
plot(results0,results1)
abline(0,1)
plot(bias0$Bias,bias1$Bias)
abline(0,1)
par(mfrow=c(1,1))

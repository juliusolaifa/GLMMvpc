glmmdata <- function(iter.size, x, beta, sigma.u, group.sizes, family,...) {
  
  x <- as.matrix(x)
  if(nrow(x) != sum(group.sizes)) {
    stop("number of rows of x must be equal to sum of cluster sizes")
  }
  X <- stats::model.matrix(~x)
  #####
  sigma.u <- as.matrix(sigma.u)
  n_groups <- length(group.sizes)
  n_rand_par <- ncol(sigma.u)
  u <- MASS::mvrnorm(n=n_groups, mu=rep(0,n_rand_par), Sigma = sigma.u)
  group_names <- paste0("group", 1:length(group.sizes))
  group_assignments <- rep(x=group_names, times=group.sizes)
  group_assignments <- factor(group_assignments)
  #####
  z0 <- stats::model.matrix(~ 0 + group_assignments)
  z <- do.call(cbind, lapply(1:ncol(u), function(i) z0 * X[, i]))
  #####
  total_sample_size <- sum(group.sizes)
  y_matrix <- matrix(nrow=iter, ncol=total_sample_size)
  #####
  #split conditional means into list according to group

  eta <- X%*%beta + z%*%c(u)
  cond_means <- split(x =exp(eta), f = rep(1:length(group.sizes), group.sizes))
  args <- list(...) 
  generate <- switch(family,
     nbinom2={
       #MASS::rnegbin
       function() unlist(mapply(MASS::rnegbin, n=group.sizes, 
                                mu=cond_means,args$theta))
     },
     tweedie={
       #tweedie::rtweedie
       function() unlist(mapply(tweedie::rtweedie, n=group.sizes, 
                                mu=cond_means, args$phi, args$power))
     },
     compois={
       #COMPoissonReg::rcmp
       function() unlist(mapply(COMPoissonReg::rcmp, n=group.sizes, 
                                lambda=cond_means, args$nu))
     },
     # genpois={
     #   #HMMpa::rgenpois
     #   lambda1 <- unlist(sapply(cond_means, function(cm) cm *(1 - args$lambda2)))
     #   function() unlist(mapply(HMMpa::rgenpois, n=group.sizes, lambda1, args$lambda2))
     # },
     stop("data generation not implemented for family: ", family)
  )
  
  for (i in 1:iter) {                
    dat <- generate()
    y_matrix[i, ] <- dat 
  }
  
  result <- structure(list("x" = x, "z" = z, "y" = y_matrix, "family" = family, 
                           "group" = group_assignments), class = "glmmData")
  result
}
# indexing method for `heritData` class
"[.glmmData" <- function(dataObj, i, ...) {
  data.frame(x = dataObj$x, group = dataObj$group, y = dataObj$y[i,])
}

print.glmmData <- function(dataObj) {
  cat("Family:", dataObj$family, "\n")
  
  nx <- ncol(dataObj$x)
  datamat <- rbind(t(dataObj$x), dataObj$y)
  colnames(datamat) <- dataObj$group
  rownames(datamat) <- c(paste0('x', 1:nx), paste0('row', 1:nrow(dataObj$y)))
  
  print(datamat)
}

summary.glmmData <- function(dataObj) {
  cat("Summary of glmmData object\n")
  cat("----------------------------\n")
  cat("Family:", dataObj$family, "\n\n")
  cat("Dimensions of data:", dim(dataObj)[1], "rows,", dim(dataObj)[2], "columns\n\n")
  cat("Number of groups:", length(unique(dataObj$group)), "\n")
  cat("Group sizes:\n")
  print(table(dataObj$group))
  invisible(dataObj)
}

dim.glmmData <- function(datObj) {
  dim(datObj$y)
}

nrow.glmmData <- function(dataObj) {
  dim(dataObj)[1]
}

head.glmmData <- function(dataObj, n = 5L, ...) {
  head_y <- utils::head(dataObj$y, n)
  structure(list(x = dataObj$x, y = head_y, family = dataObj$family, 
                 group = dataObj$group), class = "glmmData")
}

tail.glmmData <- function(dataObj, n = 5L, ...) {
  tail_y <- utils::tail(dataObj$y, n)
  structure(list(x = dataObj$x, y = tail_y, family = dataObj$family, 
                 group = dataObj$group), class = "glmmData")
}

# Example
#x <- c(23,24,25,26,29,30,18,22,21,16,18,16,21)
# x <- c(1,1,0,1,0,1,0,0,1,0,1,1,0)
# beta <- c(7,9)
# sigma.u <- matrix(c(2,1,1,2),2)
# group.sizes <- c(3,2,5,3)
# theta = 2
# phi = 2
# power = 1.6
# nu = 3
# lambda2 = 0.7
# iter =10
# data1 <- glmmdata(iter, x, beta, sigma.u, group.sizes, "nbinom2", theta=theta)
# data2 <- glmmdata(iter, x, beta, sigma.u, group.sizes, "tweedie", phi=phi, power=power)
# data3 <- glmmdata(iter, x, beta, sigma.u, group.sizes, "compois", nu=nu)
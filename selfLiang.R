library(MASS)

# Initialize matrix
I_theta <- matrix(1:16, nrow=4)
n <- dim(I_theta)[1]
mu <- rep(0,n)

# Function to compute projection matrix
computeProjectionMatrix <- function(basisMatrix) {
  I <- diag(n)
  P <- I - I_theta %*% basisMatrix %*% ginv(t(basisMatrix) %*% I_theta %*% basisMatrix) %*% t(basisMatrix)
  return(P)
}



# Initialize basis matrices
B1 <- matrix(0, nrow=n, ncol=n)
B2 <- matrix(0, nrow=n, ncol=n); B2[1,1] <- 1
B3 <- matrix(0, nrow=n, ncol=n); B3[2,2] <- 1
B4 <- matrix(0, nrow=n, ncol=n); B4[1,1] <- 1; B4[2,2] <- 1

# Compute projection matrices using the function
P1 <- computeProjectionMatrix(B1)
P2 <- computeProjectionMatrix(B2)
P3 <- computeProjectionMatrix(B3)
P4 <- computeProjectionMatrix(B4)

# Print results
print(P1 %*% I_theta %*% t(P1))
print(P2 %*% I_theta %*% t(P2))
print(P3 %*% I_theta %*% t(P3))
print(P4 %*% I_theta %*% t(P4))
print(mu)

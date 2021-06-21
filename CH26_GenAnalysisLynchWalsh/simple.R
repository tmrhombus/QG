
# cleanup if you're running interactively
rm(list=ls())

# load functions
source("usefulfunctions.R")

# MME is y=Xb+Zu+e

#load data
data <- read.csv("data.csv")

# define initial matrices
data$env <- as.factor(data$env)
data$sire <- as.factor(data$sire)

y <- data$pheno
X <- model.matrix(~env - 1, data = data)
Z <- model.matrix(~sire -1, data = data)

# Covariance Matrices, hmm
 # R = s^2_E I 6x6 environments [X]
 # G = s^2_s I 3x3 sire effects [Z]

R <- diag(6)
G <- diag(3)

sig2_E <- 6
sig2_A <- 8
sig2_S <- sig2_A/4

R <- R * sig2_E
G <- G * sig2_S


# V = covariance of y
#V = Z %*% G %*% t(Z) + R
V <- makeCOVy(Z, G, R)

# BLUE of b
bhat <- makeBLUE(X, V, y)
#bhat <- solve( t(X) %*% solve(V) %*% X ) %*% t(X) %*% solve(V) %*% y
# BLUP of u
uhat <- makeBLUP(X, Z, y, V, G)
#uhat <- G %*% t(Z) %*% solve(V) %*% ( y - X %*% bhat )

# print values
# these give the same values as in the book
# example 2
bhat
uhat

# example 4
blueblup <- makeBLUEBLUP(X, Z, y, R, G)

# checks out
blueblup



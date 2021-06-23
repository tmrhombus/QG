
rm(list=ls())

# Example 2

setwd("~/Repositories/QG/CH26_GenAnalysisLynchWalsh/")
data <- read.csv("data.csv")

data$env <- as.factor(data$env)
data$sire <- as.factor(data$sire)

y <- data$pheno
X <- model.matrix(~env - 1, data = data)
Z <- model.matrix(~sire -1, data = data)

Ve <- 6
Vs <- 8/4

G <- Vs * diag(3)
R <- Ve * diag(6)

V <- Z %*% G %*% t(Z) + R
V_inv <- solve(V)

b <- solve(t(X) %*% V_inv %*% X) %*% t(X) %*% V_inv %*% y
u <- G %*% t(Z) %*% V_inv %*% (y - X %*% b)

# > u
# [,1]
# [1,] -0.05555556
# [2,]  0.11111111
# [3,] -0.05555556


# Using MME script
source("MME-updated.R")
BLUP1 <- MME(y=y,X=X,Z=list(Z),Vu=list(G),Ve=Ve,m=1)

# > BLUP1$BLUP
# [1] 10.58333 10.75000 10.58333






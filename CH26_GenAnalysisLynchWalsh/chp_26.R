
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
mean(b) + u

# > u
# [,1]
# [1,] -0.05555556
# [2,]  0.11111111
# [3,] -0.05555556


# Using MME script
source("MME-updated.R")
BLUP1 <- MME(y=y,X=X,Z=list(Z),Vu=list(G),Ve=Ve,m=1)
BLUP1


# > BLUP1$BLUP
# [1] 10.58333 10.75000 10.58333



# MME

Q11 <- t(X) %*% solve(R) %*% X

tmp <- crossprod(X,X)
tmp/Ve
Jeff <- crossprod(X,X/Ve)

Q12 <- t(X) %*% solve(R) %*% Z
Q21 <- t(Z) %*% solve(R) %*% X
Q22 <- (t(Z) %*% solve(R) %*% Z) + solve(G)

LH <- rbind(cbind(Q11, Q12), cbind(Q21, Q22))

RH1 <- t(X) %*% solve(R) %*% y
RH2 <- t(Z) %*% solve(R) %*% y

RH <- rbind(RH1, RH2)

solve(Q) %*% RH

# [,1]
# env1   8.22222222
# env2  13.05555556
# sire1 -0.05555556
# sire2  0.11111111
# sire3 -0.05555556

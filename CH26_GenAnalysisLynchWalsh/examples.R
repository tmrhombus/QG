rm(list=ls())

# MME y=Xb+Zu+e

# define initial matrices

# y = phenotypic observations
y = matrix(
 c(9, 12, 11, 6, 7, 14),
 nrow=6,
 ncol=1,
 byrow=TRUE
)

# X = design (environment)
X = matrix(
 c(1,0, 0,1, 1,0, 1,0, 1,0, 0,1),
 nrow=6,
 ncol=2,
 byrow=TRUE
) 

# b = environments
# b = matrix(
#  c(b1, b2),
#  nrow=2,
#  ncol=1,
#  byrow=TRUE
# )

# incidence (genetic)
Z = matrix(
 c(1,0,0, 1,0,0, 0,1,0, 0,1,0, 0,0,1, 0,0,1),
 nrow=6,
 ncol=3,
 byrow=TRUE
)

# u sires
# u = matrix(
#  c(u1, u2, u3),
#  nrow=3,
#  ncol=1,
#  byrow=TRUE
# )


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
V = Z %*% G %*% t(Z) + R

# BLUE of b
bhat <- solve( t(X) %*% solve(V) %*% X ) %*% t(X) %*% solve(V) %*% y
# BLUP of u
uhat <- G %*% t(Z) %*% solve(V) %*% ( y - X %*% bhat )

# these give the same values as in the book
# example 2
bhat
uhat

# example 4
# MME equations
#
# | XT R-1 X   XT R-1 Z       | |bhat|   |XT R-1 y|
# | ZT R-1 X   ZT R-1 Z + G-1 | |uhat| = |ZT R-1 y| 
# A B = C

# left side matrix
# (top left) 
TL <- t(X) %*% solve(R) %*% X
# (top right)
TR <- t(X) %*% solve(R) %*% Z
# (bot left)
BL <- t(Z) %*% solve(R) %*% X
# (bot right)
BR <- t(Z) %*% solve(R) %*% Z + solve(G)

# # check
#  6 * TL
#  6 * TR
#  6 * BL
#  6 * BR

# right side matrix
RT <- t(X) %*% solve(R) %*% y
RB <- t(Z) %*% solve(R) %*% y

# # check
#  6 * RT
#  6 * RB
 
# this doesn't work
# # MME : LHS * bhat/uhat = RHS
# LHS <- matrix(
#  c(TL, TR, BL, BR),
#   nrow=5,
#   ncol=5,
#   byrow=TRUE
# )

# binding a matrix from submatrices
LHS_T <- cbind(TL,TR)
LHS_B <- cbind(BL,BR)
LHS <- rbind(LHS_T,LHS_B)

RHS <- rbind(RT,RB)

# B = A-1 C
bhatuhat <- solve(LHS) %*% RHS

# checks out
bhatuhat



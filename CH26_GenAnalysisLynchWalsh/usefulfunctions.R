
makeCOVy <- function(Z, G, R){
 # Z 
 V <-  Z %*% G %*% t(Z) + R
 return (V)
}
 
makeBLUE <- function(X, V, y){
 bhat <- solve( t(X) %*% solve(V) %*% X ) %*% t(X) %*% solve(V) %*% y
 return (bhat)
}
 
 
makeBLUP <- function(X, Z, y, V, G){
 uhat <- G %*% t(Z) %*% solve(V) %*% ( y - X %*% makeBLUE(X, V, y) )
 return (uhat)
}


makeBLUEBLUP <- function(X, Z, y, R, G){
 # MME equations
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
 
 # RIGHT side matrix
 RT <- t(X) %*% solve(R) %*% y
 RB <- t(Z) %*% solve(R) %*% y
 RHS <- rbind(RT,RB)
 
 # binding LEFT matrix from submatrices
 LHS_T <- cbind(TL,TR)
 LHS_B <- cbind(BL,BR)
 LHS <- rbind(LHS_T,LHS_B)
 
 # B = A-1 C
 blueblup <- solve(LHS) %*% RHS
 
 return <- (blueblup)
}

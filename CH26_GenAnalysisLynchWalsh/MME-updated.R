direct.sum <- function(argv) {
  i = 0
  for( a in argv ) {
    m = as.matrix(a)
    if(i == 0)
      rmat = m
    else
    {
      nr = dim(m)[1]
      nc = dim(m)[2]
      aa = cbind(matrix(0,nr,dim(rmat)[2]),m)
      rmat = cbind(rmat,matrix(0,dim(rmat)[1],nc))
      rmat = rbind(rmat,aa)
    }
    i = i+1
  }
  return(rmat)
}

MME <- function(y,X,Z,Vu,Ve,m=1) {
  #y = vector of phenotypes (without missing data)
  #X = design matrix for fixed effects, assumes single factor.
  #Z = list of design matrices for random effects. Function returns BLUP for first random effect.
  #Vu = list of var-covar matrices for the random effects.
  #Ve = vector of residual variances, same length as y
  #Returns BLUP and r2 for the sum of first m random effects for each genotype
  stopifnot(all(!is.na(y)))
  gid <- gsub("Name","",colnames(Z[[1]]))
  n <- length(gid)
  Z2 <- NULL
  for (i in 1:length(Z)) {
    Z2 <- cbind(Z2,Z[[i]])
  }
  Z <- Z2
  Q11 <- crossprod(X,X/Ve)
  Q12 <- crossprod(X,Z/Ve)
  Q21 <- t(Q12)
  Q22 <- crossprod(Z,Z/Ve)+direct.sum(lapply(Vu,solve))
  Q <- rbind(cbind(Q11,Q12),cbind(Q21,Q22))
  b1 <- crossprod(X,y/Ve)
  b2 <- crossprod(Z,y/Ve)
  b <- c(b1,b2)
  C <- solve(Q)
  ans <- C%*%b
  n.env <- ncol(X)
  BLUP <- mean(ans[1:n.env]) #fixed effect
  for (i in 1:m) {
    BLUP <- BLUP + ans[n.env+n*(i-1)+1:n]
  }
  ix <- n.env+1:(n*m)
  PEV <- C[ix,ix]
  v <- NULL
  denom <- 0
  for (i in 1:m) {
    v <- cbind(v,diag(n))
    denom <- denom + diag(Vu[[i]])
  }
  r2 <- 1 - diag(v%*%PEV%*%t(v))/denom
  tmp <- data.frame(Name=gid,BLUP=BLUP,r2=r2,stringsAsFactors = F)
  rownames(tmp) <- NULL
  return(tmp)
}



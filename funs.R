library(Matrix)
library(glmgen)

Seq = function(a,b) {
  if (a<=b) return(a:b)
  else return(integer(0))
}

# Construct 1d difference operator
get.1d.sparse = function(n) {
  return(bandSparse(n-1, m=n, k=c(0, 1), 
                    diagonals=list(rep(-1,n-1),rep(1,n-1))))
}

# Construct 2d difference operator
get.2d.sparse = function(n1, n2) {
  D1 = bandSparse(n1 * n2, m = n1 * n2, k = c(0, 1), 
    diagonals = list(rep(-1, n1 * n2), rep(1, n1 * 
      n2 - 1)))
  D1 = D1[(Seq(1, n1 * n2)%%n1) != 0, ]
  D2 = bandSparse(n1 * n2 - n1, m = n1 * n2,
    k = c(0, n1), diagonals = list(rep(-1, n1 * n2),
                      rep(1, n1 * n2 - 1)))
  return(rBind(D1, D2))
}

crit = function(y, b, D, lambda) {
  return(1/2*sum((y-b)^2) + lambda*sum(abs(D%*%b)))
}

# Box function
bx = function(x, s) {
  y = x
  y[x > s] = s
  y[x < -s] = -s
  return(y)
}

# Proximal gradient descent on dual
prox.grad.genlasso = function(y, D, lambda,
  u0=rep(0,nrow(D)), step=1, maxsteps=2000,
  verbose=FALSE) {

  y = as.numeric(y)
  u = u0
  for (k in 1:maxsteps) {
    if (verbose && k%%100==0) cat(k,"...")
    u = bx(as.numeric(u+step*D%*%(y-t(D)%*%u)),lambda)
  }
  if (verbose) cat("\n")
  return(list(u=u,b=as.numeric(y-t(D)%*%u)))
}

# Coordinate descent on dual
coord.desc.genlasso = function(y, D, lambda,
  u0=rep(0,nrow(D)), maxsteps=100, verbose=FALSE) {

  y = as.numeric(y)
  m = nrow(D)
  u = u0
  DDt = crossprod(t(D))
  Dy = as.numeric(D%*%y)
  for (k in 1:maxsteps) {
    if (verbose && k%%10==0) cat(k,"...")
    a = DDt %*% u
    for (i in 1:m) {
      if (i > 1) a = aold - DDt[,i]*uold + DDt[,i]*u[i]
      aold = a
      uold = u[i]
      u[i] = bx((Dy[i] - a[i])/DDt[i,i], lambda)
    }
  }
  if (verbose) cat("\n")
  return(list(u=u,b=as.numeric(y-t(D)%*%u)))
}

# Ryan's implementation of specialized ADMM for the 2d fused lasso
special.admm.2d = function(ymat, lambda, rho=1, maxsteps=25,
  verbose=FALSE) {

  n1 = nrow(ymat)
  n2 = ncol(ymat)
  x = z = u = matrix(0,n1,n2)

  for (k in 1:maxsteps) {
    if (verbose) cat(k,"... ")
    # x-update
    for (i in 1:n1) {
      x[i,] = trendfilter((ymat[i,]+rho*(z[i,]-u[i,]))/(1+rho),
         k=0, lambda=2*lambda/(1+rho))$beta
    }
    # z-update
    for (j in 1:n2) {
      z[,j] = trendfilter((ymat[,j]+rho*(x[,j]+u[,j]))/(1+rho),
         k=0, lambda=2*lambda/(1+rho))$beta
    }
    # u-update
    u = u + x - z
  }

  return(x)
}

# Chambolle and Darbon's "exact" parametric max flow algortihm
tv.exact = function(y, lambda, err=0) {
  sy = ncol(y)
  sx = nrow(y)
  
  dyn.load("tvh-exact.so")

  o = .C("TV4_R",
    y=as.double(y),
    sx=as.integer(sx),
    sy=as.integer(sy),
    lambda=as.double(lambda),
    err=as.double(err),
    dup=FALSE)

  return(matrix(o$y,sx,sy))
}

# Chain approximation to the 2d fused lasso problem
chain.approx.2d = function(ymat, lambda, niter=0, verbose=FALSE) {
  n1 = nrow(ymat); n2 = ncol(ymat)
  
  # Chain by columns
  b = trendfilter(chain.mat.cols(ymat), k=0, lambda=lambda)$beta
  beta.cols = unchain.mat.cols(b,n1,n2)
  
  # Chain by rows
  b = trendfilter(chain.mat.rows(ymat), k=0, lambda=lambda)$beta
  beta.rows = unchain.mat.rows(b,n1,n2)
  
  # Average the two
  beta.avg = (beta.cols+beta.rows)/2

  # Iterative
  beta.iter = NULL
  if (niter > 0) {
    beta.iter = array(dim=c(n1,n2,niter))
    beta.iter[,,1] = beta.cols
    rows = TRUE
    for (k in 2:niter) {
      if (verbose) cat(k,"... ")
      if (rows) {
        b = trendfilter(chain.mat.rows(beta.iter[,,k-1]), k=0, lambda=lambda)$beta
        beta.iter[,,k] = unchain.mat.rows(b,n1,n2)
      }
      else {
        b = trendfilter(chain.mat.cols(beta.iter[,,k-1]), k=0, lambda=lambda)$beta
        beta.iter[,,k] = unchain.mat.cols(b,n1,n2)
      }
      rows = !rows
    }
  }
  
  return(list(beta.avg=beta.avg, beta.cols=beta.cols, beta.rows=beta.rows,
              beta.iter=beta.iter))
}

chain.approx.2d.cd = function(ymat, lambda, niter=10, verbose=FALSE) {
  n1 = nrow(ymat); n2 = ncol(ymat)
  y.cols = chain.mat.cols(ymat)
  y.rows = chain.mat.rows(ymat)
  beta.iter = array(dim=c(n1,n2,niter))
  
  D = get.1d.sparse(n1*n2)
  DDt = crossprod(t(D))

  # Chain first by columns
  b = trendfilter(chain.mat.cols(ymat), k=0, lambda=lambda)$beta
  u = solve(DDt, D%*%(y.cols-b))
  beta.iter[,,1] = unchain.mat.cols(b,n1,n2)

  rows = TRUE
  for (k in Seq(2,niter)) {
    if (verbose) cat(k,"... ")
    if (rows) {
      # Note that b solved over columns
      # Set to 0 all entries of u corresponding to shared edges
      # So every (n1)st entry
      u[seq(n1,length(u),by=n1)] = 0

      # Define residual, solve for b over rows
      r = chain.mat.rows(unchain.mat.cols(as.numeric(y.cols-t(D)%*%u),n1,n2))
      b = trendfilter(r, k=0, lambda=lambda)$beta
      u = solve(DDt, D%*%(r-b))
      beta.iter[,,k] = unchain.mat.rows(b,n1,n2)
    }
    else {
      # Note that b solved over rows
      # Set to 0 all entries of u corresponding to shared edges
      # So every (n2)st entry
      u[seq(n2,length(u),by=n2)] = 0

      # Define residual, solve for b over rows
      r = chain.mat.cols(unchain.mat.rows(as.numeric(y.rows-t(D)%*%u),n1,n2))
      b = trendfilter(r, k=0, lambda=lambda)$beta
      u = solve(DDt, D%*%(r-b))
      beta.iter[,,k] = unchain.mat.cols(b,n1,n2)
    }
    rows = !rows
  }

  return(list(beta=beta.iter[,,niter], beta.iter=beta.iter))
}

chain.mat.cols = function(a) {
  even.cols = seq(2, ncol(a), by=2)
  a[,even.cols] = apply(a[,even.cols,drop=F], 2, rev)
  return(as.numeric(a))
}

unchain.mat.cols = function(a, n1, n2) {
  a = matrix(a, n1, n2)
  even.cols = seq(2, n2, by=2)
  a[,even.cols] = apply(a[,even.cols,drop=F], 2, rev)
  return(a)
}

chain.mat.rows = function(a) {
  even.rows = seq(2, nrow(a), by=2)
  a[even.rows,] = t(apply(a[even.rows,,drop=F], 1, rev))
  return(as.numeric(t(a)))
}

unchain.mat.rows = function(a, n1, n2) {
  a = matrix(a, n1, n2, byrow=TRUE)
  even.rows = seq(2, n1, by=2)
  a[even.rows,] = t(apply(a[even.rows,,drop=F], 1, rev))
  return(a)
}

# Make "C" image
make.c = function(n1, n2, vals=1:4) {
  k1 = n1/2
  k2 = n2/2
  y = matrix(vals[1],2*k1,2*k2)
  y[1:k1,1:k2] = y[k1:(2*k1),k2:(2*k2)] = vals[1]
  y[1:k1,k2:(2*k2)] = y[k1:(2*k1),1:k2] = vals[2]
  # Circles
  for (i in (0.5*k1):(1.5*k1)) {
    for (j in 1:k2) {
      if (sqrt((i-k1)^2+(j-0.5*k2)^2)<0.5*min(k1,k2)) {
        y[i,j] = vals[3]
      }
    }
  }
  for (i in k1:(2*k1)) {
    for (j in k2:(2*k2)) {
      if (sqrt((i-2*k1)^2+(j-2*k2)^2)<0.9*min(k1,k2)) {
        y[i,j] = vals[3]
      }
    }
  } 
  # CMU C
  y[(0.7*k1):(1.3*k1),(1.5*k2):(1.75*k2)] = vals[4]
  y[(0.7*k1):(1.3*k1),(0.5*k2):(0.75*k2)] = vals[4]
  y[(0.5*k1):(0.7*k1),(0.75*k2):(1.5*k2)] = vals[4]
  y[(1.3*k1):(1.5*k1),(0.75*k2):(0.9*k2)] = vals[4]
  y[(1.3*k1):(1.5*k1),(1.35*k2):(1.5*k2)] = vals[4]
  x = upper.tri(matrix(nrow=length((1.3*k1):(1.5*k1)),ncol=length((0.5*k2):(0.75*k2))),diag=TRUE)
  y[(1.3*k1):(1.5*k1),(0.5*k2):(0.75*k2)][x] = vals[4]
  x = upper.tri(matrix(nrow=length((0.5*k1):(0.7*k1)),ncol=length((0.5*k2):(0.75*k2))),diag=TRUE)
  y[(0.5*k1):(0.7*k1),(0.5*k2):(0.75*k2)][x[nrow(x):1,]] = vals[4]
  x = upper.tri(matrix(nrow=length((1.3*k1):(1.5*k1)),ncol=length((1.5*k2):(1.75*k2))),diag=TRUE)
  y[(1.3*k1):(1.5*k1),(1.5*k2):(1.75*k2)][x[,ncol(x):1]] = vals[4]
  x = upper.tri(matrix(nrow=length((0.5*k1):(0.7*k1)),ncol=length((1.5*k2):(1.75*k2))),diag=TRUE)
  y[(0.5*k1):(0.7*k1),(1.5*k2):(1.75*k2)][x[nrow(x):1,ncol(x):1]] = vals[4]
  x = upper.tri(matrix(nrow=length((0.5*k1):(0.7*k1)),ncol=length(k2:(1.25*k2))),diag=TRUE)
  y[(0.5*k1):(0.7*k1),k2:(1.25*k2)][x[nrow(x):1,]] = vals[4]
  return(y)
}


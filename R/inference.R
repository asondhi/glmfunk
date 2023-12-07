# Debias test for Poisson regression with individual intercepts alpha; theta = (alpha, beta)'
debias.poisson.funk = function(y, X, thetahat, mu) {
  
  n = length(y)
  p = ncol(X)
  Xtilde = cbind(diag(n), X)
  betahat = thetahat[-(1:n)]
  
  # Gradient term
  mu.est = exp(as.numeric(Xtilde%*%thetahat))
  grad.est = as.numeric(t(X)%*%(y - mu.est))/n
  
  # Hessian term
  W.est = (y - mu.est) %*% t(y - mu.est) # sandwich
  sigmahat.est = t(X)%*%W.est%*%X/n
  
  M.est = InverseLinfty(sigma = sigmahat.est, n = n, resol=1.5, mu=mu, maxiter=50, threshold=1e-2, verbose = FALSE)

  # Asymptotic covariance matrix 
  cov.est = M.est%*%sigmahat.est%*%t(M.est)
  
  beta.db = betahat + M.est%*%grad.est
  ts = as.numeric(sqrt(n)*beta.db/sqrt(diag(cov.est)))
  return(ts)
}

debias.poisson.funk.ci = function(y, X, thetahat, mu) {
  
  n = length(y)
  p = ncol(X)
  Xtilde = cbind(diag(n), X)
  betahat = thetahat[-(1:n)]
  
  # Gradient term
  mu.est = exp(as.numeric(Xtilde%*%thetahat))
  grad.est = as.numeric(t(X)%*%(y - mu.est))/n
  
  # Hessian term
  W.est = (y - mu.est) %*% t(y - mu.est) # sandwich
  sigmahat.est = t(X)%*%W.est%*%X/n
  
  M.est = InverseLinfty(sigma = sigmahat.est, n = n, resol=1.5, mu=mu, maxiter=50, threshold=1e-2, verbose = FALSE)
  
  # Asymptotic covariance matrix 
  cov.est = M.est%*%sigmahat.est%*%t(M.est)
  
  beta.db = betahat + M.est%*%grad.est
  
  ci_lower = beta.db - 1.96 * sqrt(diag(cov.est))/sqrt(n)
  ci_upper = beta.db + 1.96 * sqrt(diag(cov.est))/sqrt(n)
  ci = as.data.frame(cbind(ci_lower, ci_upper))
  
  return(ci)
}


# Debias test for Poisson regression with ordinary lasso
debias.poisson.lasso = function(y, X, betahat) {
  
  n = length(y)
  p = ncol(X)

  # Gradient term
  mu.lasso = exp(as.numeric(X%*%betahat))
  grad.lasso = as.numeric(t(X)%*%(y - mu.lasso))/n

  # Hessian term
  #  W.lasso = diag(mu.lasso)
  W.lasso = (y - mu.lasso) %*% t(y - mu.lasso) # sandwich
  sigmahat.lasso = t(X)%*%W.lasso%*%X/n
  M.lasso = InverseLinfty(sigma = sigmahat.lasso, n = n, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = FALSE)

  # Asymptotic covariance matrix 
  cov.est = M.lasso%*%sigmahat.lasso%*%t(M.lasso)

  beta.db = betahat + M.lasso%*%grad.lasso
  ts = as.numeric(sqrt(n)*beta.db/sqrt(diag(cov.est)))
  return(ts)
}

debias.poisson.lasso.ci = function(y, X, betahat) {
  
  n = length(y)
  p = ncol(X)
  
  # Gradient term
  mu.lasso = exp(as.numeric(X%*%betahat))
  grad.lasso = as.numeric(t(X)%*%(y - mu.lasso))/n
  
  # Hessian term
  #  W.lasso = diag(mu.lasso)
  W.lasso = (y - mu.lasso) %*% t(y - mu.lasso) # sandwich
  sigmahat.lasso = t(X)%*%W.lasso%*%X/n
  M.lasso = InverseLinfty(sigma = sigmahat.lasso, n = n, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = FALSE)
  
  # Asymptotic covariance matrix 
  cov.est = M.lasso%*%sigmahat.lasso%*%t(M.lasso)
  
  beta.db = betahat + M.lasso%*%grad.lasso
  
  ci_lower = beta.db - 1.96 * sqrt(diag(cov.est))/sqrt(n)
  ci_upper = beta.db + 1.96 * sqrt(diag(cov.est))/sqrt(n)
  ci = as.data.frame(cbind(ci_lower, ci_upper))
  
  return(ci)
}


# Debias test for logistic regression with individual intercepts alpha; theta = (alpha, beta)'
debias.logistic.funk = function(y, X, thetahat) {
  
  n = length(y)
  p = ncol(X)
  Xtilde = cbind(diag(n), X)
  betahat = thetahat[-(1:n)]

  # Gradient term
  p.est = plogis(as.numeric(Xtilde%*%thetahat))
  grad.est = as.numeric(t(X)%*%(y - p.est))/n

  # Hessian term
  W.est = diag(p.est*(1-p.est))
  sigmahat.est = t(X)%*%W.est%*%X/n
  M.est = InverseLinfty(sigma = sigmahat.est, n = n, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = FALSE)

  # Asymptotic covariance matrix 
  cov.est = M.est%*%sigmahat.est%*%t(M.est)

  beta.db = betahat + M.est%*%grad.est
  ts = as.numeric(sqrt(n)*beta.db/sqrt(diag(cov.est)))
  return(ts)
}

# Debias test for logistic regression with ordinary lasso
debias.logistic.lasso = function(y, X, betahat) {
  
  n = length(y)
  p = ncol(X)

  # Gradient term
  p.lasso = plogis(as.numeric(X%*%betahat))
  grad.lasso = as.numeric(t(X)%*%(y - p.lasso))/n

  # Hessian term
  W.lasso = diag(p.lasso*(1-p.lasso))
  sigmahat.lasso = t(X)%*%W.lasso%*%X/n
  M.lasso = InverseLinfty(sigma = sigmahat.lasso, n = n, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = FALSE)

  # Asymptotic covariance matrix 
  cov.est = M.lasso%*%sigmahat.lasso%*%t(M.lasso)

  beta.db = betahat + M.lasso%*%grad.lasso
  ts = as.numeric(sqrt(n)*beta.db/sqrt(diag(cov.est)))
  return(ts)
}

# Debias test for linear regression with individual intercepts alpha; theta = (alpha, beta)'
debias.linear.funk = function(y, X, thetahat) {
  
  n = length(y)
  p = ncol(X)
  Xtilde = cbind(diag(n), X)
  betahat = thetahat[-(1:n)]

  # Gradient term
  grad.est = as.numeric(t(X)%*%(y - as.numeric(Xtilde%*%thetahat)))/n

  # Hessian term
  sigmahat.est = t(X)%*%X/n
  M.est = InverseLinfty(sigma = sigmahat.est, n = n, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = FALSE)

  # Asymptotic covariance matrix 
  cov.est = M.est%*%sigmahat.est%*%t(M.est)

  beta.db = betahat + M.est%*%grad.est
  
  # Estimate Gaussian noise parameter
#  scal_m = scalreg(X, y)
  scal_m = scalreg(Xtilde, y)
  s.hat = scal_m$hsigma
  
  ts = as.numeric(sqrt(n)*beta.db/(s.hat * sqrt(diag(cov.est))))
  return(ts)
}

# Debias test for linear regression with ordinary lasso
debias.linear.lasso = function(y, X, betahat) {
  
  n = length(y)
  p = ncol(X)

  # Gradient term
  grad.lasso = as.numeric(t(X)%*%(y - as.numeric(X%*%betahat)))/n

  # Hessian term
  sigmahat.lasso = t(X)%*%X/n
  M.lasso = InverseLinfty(sigma = sigmahat.lasso, n = n, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = FALSE)

  # Asymptotic covariance matrix 
  cov.est = M.lasso%*%sigmahat.lasso%*%t(M.lasso)

  beta.db = betahat + M.lasso%*%grad.lasso
  
  # Estimate Gaussian noise parameter
  scal_m = scalreg(X, y)
  s.hat = scal_m$hsigma
  
  ts = as.numeric(sqrt(n)*beta.db/(s.hat * sqrt(diag(cov.est))))
  return(ts)
}


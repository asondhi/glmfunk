# Functions to fit funk-l1 models: training, prediction, cross-validation

funkl1.fit = function(y, X, Ln, Lp, reg_params, model = c("linear", "logistic", "poisson"), mu = 1e-03, theta_init = NULL, tol = 1e-06, max_iter = 1e6, verbose = FALSE) {
  
  n = nrow(X)
  p = ncol(X)
  y_m = as.matrix(y, ncol = 1)
  if (is.null(theta_init)) {
    theta_init = rep(0, n + p)
  }
  theta_init_m = as.matrix(theta_init, ncol = 1)
  
  gamma_n = reg_params[1] # subject graph regularization of alphas 
  gamma_p = reg_params[2] # feature graph regularization of betas
  lambda = reg_params[3] # lasso regularization of betas 
  Ln = Ln + 0.01*diag(n)
  
  Xtilde = cbind(diag(n), X)
  
  if (model == "poisson") {
    XX_eigen = max(eigen(t(Xtilde)%*%Xtilde)$values)
    Lpart = base::norm(gamma_n*Ln, "2") + base::norm(gamma_p*Lp, "2")^2/mu
    thetahat = as.numeric(funk_l1_poisson_fit(y_m, X, gamma_n*Ln, gamma_p*Lp, lambda, mu, XX_eigen, Lpart, theta_init_m, tol, max_iter, verbose))
  }
  
  if (model == "logistic") {
    L = 0.25*max(eigen(t(Xtilde)%*%Xtilde)$values) + base::norm(gamma_n*Ln, "2") + base::norm(gamma_p*Lp, "2")^2/mu
    thetahat = as.numeric(funk_l1_logistic_fit(y_m, X, gamma_n*Ln, gamma_p*Lp, lambda, mu, L, theta_init_m, tol, max_iter, verbose))
  }
  
  if (model == "linear") {
    L = max(eigen(t(Xtilde)%*%Xtilde)$values) + base::norm(gamma_n*Ln, "2") + base::norm(gamma_p*Lp, "2")^2/mu
    thetahat = as.numeric(funk_l1_linear_fit(y_m, X, gamma_n*Ln, gamma_p*Lp, lambda, mu, L, theta_init_m, tol, max_iter, verbose))
  }  
  
  alphahat = thetahat[1:n]
  betahat = thetahat[(n+1):(n+p)]
  
  if (any(is.na(alphahat))) {
    alphahat = rep(0, n)
    warning("did not converge, giving alphahat as zeroes")
  }
  if (any(is.na(betahat))) {
    betahat = rep(0, p)
    warning("did not converge, giving betahat as zeroes")
  }
  
  output = list(betahat = betahat, alphahat = alphahat, model = model)
  return(output)
}

funkl1.predict = function(funkout, Xtest, Ln_full, train_ind) {
  
  nfull = nrow(Ln_full)
  test_ind = setdiff(1:nfull, train_ind)
  
  # Predict node effects for test subjects 
  Ln_22 = Ln_full[test_ind, test_ind]
  Ln_21 = Ln_full[test_ind, train_ind]
  alphatrain = as.matrix(funkout$alphahat, ncol=1)
  alphatest = as.numeric(-ginv(Ln_22) %*% Ln_21 %*% alphatrain)
  
  eta.test = alphatest + as.numeric(Xtest%*%funkout$betahat)
  
  if (funkout$model == "poisson") {
    out = exp(eta.test)
  }  
  
  if (funkout$model == "logistic") {
    out = plogis(eta.test)
  }
  
  if (funkout$model == "linear") {
    out = eta.test
  }
  
  return(out)
}

poisson_loss = function(y, yhat) {
  mean( yhat - y * log(yhat) )
} 

logistic_dev = function(y, probs) {
  probs[probs > (1-1e-6)] = 1-(1e-6)
  probs[probs < 1e-6] = 1e-6
  mean( y*log(1/probs) + (1-y)*log(1/(1 - probs)) )
}

# CV for a single training-testing split over a parameter grid, with warm starts
# Grid should have two parameters constant, third one varying and descending 
funkl1.cv.single = function(y, X, Ln, Lp, reg_grid, train_ind, model = c("linear", "logistic", "poisson")) {
  
  n = nrow(X)
  test_ind = setdiff(1:n, train_ind)
  ytrain = y[train_ind]
  Xtrain = X[train_ind,]
  Lntrain = Ln[train_ind, train_ind]
  ytest = y[test_ind]
  Xtest = X[test_ind,]
  
  err.ests = numeric(nrow(reg_grid))
  for (i in 1:nrow(reg_grid)) {
    
    reg_i = reg_grid[i,]
    if (i == 1) {
      funk.m = funkl1.fit(ytrain, Xtrain, Lntrain, Lp, reg_i, model)
      theta_ws = c(funk.m$alphahat, funk.m$betahat)
    }
    else {
      funk.m = funkl1.fit(ytrain, Xtrain, Lntrain, Lp, reg_i, model, theta_init = theta_ws)
    }
    pred = funkl1.predict(funk.m, Xtest, Ln, train_ind)
    
    if (model == "poisson") {
      err.ests[i] = poisson_loss(ytest, pred)
    }
    
    if (model == "logistic") {
      err.ests[i] = logistic_dev(ytest, pred)
    }
    
    if (model == "linear") {
      err.ests[i] = mean((ytest - pred)^2)
    }
    
  }
  return(err.ests)  
}

# Wrap around above, do K times
funkl1.cv.Kfold = function(y, X, Ln, Lp, reg_grid, foldind, model = c("linear", "logistic", "poisson")) {
  n = nrow(X)
  nfolds = max(foldind)
  cv.ests = matrix(nrow = nrow(reg_grid), ncol = nfolds)
  for (k in 1:nfolds) {
    train_ind = which(foldind != k)
    cv.ests[,k] = funkl1.cv.single(y, X, Ln, Lp, reg_grid, train_ind, model)
  }
  return(rowMeans(cv.ests))
}

funkl1.cv.coorddesc = function(y, X, Ln, Lp, param_list, model = c("linear", "logistic", "poisson"), nfolds = 5, verb.cv = FALSE) {
  # coordinate descent with warm starts
  n = nrow(X)
  foldind = sample(nfolds, size = n, replace = TRUE) # change to be exactly n/K equal?
  
  gamma_n_grid = sort(param_list[[1]], decreasing = TRUE)
  gamma_p_grid = sort(param_list[[2]], decreasing = TRUE)
  lambda_grid = sort(param_list[[3]], decreasing = TRUE)
  
  gamma_n = gamma_n_grid[1]
  gamma_p = gamma_p_grid[1]
  lambda = lambda_grid[1]
  reg_params_curr = c(gamma_n, gamma_p, lambda)
  cvg = FALSE
  iter = 0
  while (!cvg) {
    # Optimize gamma_n
    reg_grid = cbind(gamma_n_grid, gamma_p, lambda)
    cvKerr = funkl1.cv.Kfold(y, X, Ln, Lp, reg_grid, foldind, model)
    gamma_n = gamma_n_grid[which.min(cvKerr)]
    if (verb.cv) print("# Optimized gamma_n #")
    # Optimize gamma_p
    reg_grid = cbind(gamma_n, gamma_p_grid, lambda)
    cvKerr = funkl1.cv.Kfold(y, X, Ln, Lp, reg_grid, foldind, model)
    gamma_p = gamma_p_grid[which.min(cvKerr)]
    if (verb.cv) print("# Optimized gamma_p #")
    # Optimize lambda
    reg_grid = cbind(gamma_n, gamma_p, lambda_grid)
    cvKerr = funkl1.cv.Kfold(y, X, Ln, Lp, reg_grid, foldind, model)
    lambda = lambda_grid[which.min(cvKerr)]
    if (verb.cv) print("# Optimized lambda #")
    # Check convergence
    iter = iter + 1
    reg_params_new = c(gamma_n, gamma_p, lambda)
    step = sum(abs(reg_params_new-reg_params_curr))
    if (verb.cv) {
      print(paste0("### Finished coordinate descent iter ", iter, " ###"))
      print(paste0("Old params: ", paste0(reg_params_curr, collapse = " ")))
      print(paste0("New params: ", paste0(reg_params_new, collapse = " ")))
    }
    reg_params_curr = reg_params_new
    cvg = (step < 1e-04 || iter >= 10)
  }
  minerr = funkl1.cv.Kfold(y, X, Ln, Lp, t(as.matrix(reg_params_curr)), foldind, model)
  return(list(reg_params_curr, minerr))
}

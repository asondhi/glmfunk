# Functions to fit genridge models: training, prediction, cross-validation, inference via debiasing

gracel2.fit = function(y, X, Lp, reg_params, model = c("linear", "logistic"), beta_init = NULL, tol = 1e-08, max_iter = 1e5, verbose = FALSE) {
  
  n = nrow(X)
  p = ncol(X)
  y_m = as.matrix(y, ncol = 1)
  if (is.null(beta_init)) {
    beta_init = runif(p)
  }
  beta_init_m = as.matrix(beta_init, ncol = 1)
  
  gamma_p = reg_params[1] # feature graph regularization of betas
  lambda = reg_params[2] # lasso regularization of betas 
  
  
  if (model == "logistic") {
    betahat = as.numeric(grace_l2_logistic_fit(y_m, X, gamma_p*Lp, lambda, beta_init_m, tol, max_iter, verbose))
  }
  if (model == "linear") {
    betahat = as.numeric(grace_l2_linear_fit(y_m, X, gamma_p*Lp, lambda, beta_init_m, tol, max_iter, verbose))
  }
  
  if (any(is.na(betahat))) {
    betahat = rep(0, p)
    warning("did not converge, giving betahat as zeroes")
  }
  
  out = list(betahat = betahat, model = model)
  return(out)
}

gracel2.predict = function(funkout, Xtest) {
  
  eta.test = as.numeric(Xtest%*%funkout$betahat)
  if (funkout$model == "logistic") {
    out = plogis(eta.test)
  }
  if (funkout$model == "linear") {
    out = eta.test
  }
  
  return(out)
}

logistic_dev = function(y, probs) {
  probs[probs > (1-1e-6)] = 1-(1e-6)
  probs[probs < 1e-6] = 1e-6
  mean( y*log(1/probs) + (1-y)*log(1/(1 - probs)) )
}

# CV for a single training-testing split over a parameter grid, with warm starts
# Grid should have one parameter constant, second one varying and descending
gracel2.cv.single = function(y, X, Lp, reg_grid, train_ind, model = c("linear", "logistic")) {
  
  n = nrow(X)
  test_ind = setdiff(1:n, train_ind)
  ytrain = y[train_ind]
  Xtrain = X[train_ind,]
  ytest = y[test_ind]
  Xtest = X[test_ind,]
  
  err.ests = numeric(nrow(reg_grid))
  for (i in 1:nrow(reg_grid)) {
    reg_i = reg_grid[i,]
    if (i == 1) {
      grace.m = gracel2.fit(ytrain, Xtrain, Lp, reg_i, model)
      beta_ws = grace.m$betahat
    }
    else {
      grace.m = gracel2.fit(ytrain, Xtrain, Lp, reg_i, model, beta_init = beta_ws)
    }
    pred = gracel2.predict(grace.m, Xtest)

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
gracel2.cv.Kfold = function(y, X, Lp, reg_grid, foldind, model = c("linear", "logistic")) {
  n = nrow(X)
  nfolds = max(foldind)
  cv.ests = matrix(nrow = nrow(reg_grid), ncol = nfolds)
  for (k in 1:nfolds) {
    train_ind = which(foldind != k)
    cv.ests[,k] = gracel2.cv.single(y, X, Lp, reg_grid, train_ind, model)
  }
  return(rowMeans(cv.ests))
}

gracel2.cv.coorddesc = function(y, X, Lp, param_list, model = c("linear", "logistic"), nfolds = 5, verb.cv = FALSE) {
  # coordinate descent with warm starts
  n = nrow(X)
  foldind = sample(nfolds, size = n, replace = TRUE) # change to be exactly n/K equal?
  
  gamma_p_grid = sort(param_list[[1]], decreasing = TRUE)
  lambda_grid = sort(param_list[[2]], decreasing = TRUE)
  
  gamma_p = gamma_p_grid[1]
  lambda = lambda_grid[1]
  reg_params_curr = c(gamma_p, lambda)
  cvg = FALSE
  iter = 0
  while (!cvg) {
    # Optimize gamma_p
    reg_grid = cbind(gamma_p_grid, lambda)
    cvKerr = gracel2.cv.Kfold(y, X, Lp, reg_grid, foldind, model)
    gamma_p = gamma_p_grid[which.min(cvKerr)]
    if (verb.cv) print("# Optimized gamma_p #")
    # Optimize lambda
    reg_grid = cbind(gamma_p, lambda_grid)
    cvKerr = gracel2.cv.Kfold(y, X, Lp, reg_grid, foldind, model)
    lambda = lambda_grid[which.min(cvKerr)]
    if (verb.cv) print("# Optimized lambda #")
    # Check convergence
    iter = iter + 1
    reg_params_new = c(gamma_p, lambda)
    step = sum(abs(reg_params_new-reg_params_curr))
    if (verb.cv) {
      print(paste0("### Finished coordinate descent iter ", iter, " ###"))
      print(paste0("Old params: ", paste0(reg_params_curr, collapse = " ")))
      print(paste0("New params: ", paste0(reg_params_new, collapse = " ")))
    }
    reg_params_curr = reg_params_new
    cvg = (step < 1e-04 || iter >= 10)
  }
  minerr = gracel2.cv.Kfold(y, X, Lp, t(as.matrix(reg_params_curr)), foldind, model)
  return(list(reg_params_curr, minerr))
}

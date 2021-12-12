# Functions to fit rnc-lasso models: training, prediction, cross-validation, inference via debiasing

rnc.lasso.fit = function(y, X, Ln, reg_params, model = c("linear", "logistic", "poisson"), theta_init = NULL, tol = 1e-08, max_iter = 1e7, verbose = FALSE) {
  
  n = nrow(X)
  p = ncol(X)
  y_m = as.matrix(y, ncol = 1)
  if (is.null(theta_init)) {
    theta_init = runif(n + p)
  }
  theta_init_m = as.matrix(theta_init, ncol = 1)
  
  gamma_n = reg_params[1] # subject graph regularization of alphas 
  lambda = reg_params[2] # lasso regularization of betas 
  Ln = Ln + 0.01*diag(n)
  
  if (model == "poisson") {
    Xtilde = cbind(diag(n), X)
    XX_eigen = max(eigen(t(Xtilde)%*%Xtilde)$values)
    
    Lpart = base::norm(Ln, "2")
    
    thetahat = as.numeric(rnc_lasso_poisson_fit(y_m, X, gamma_n*Ln, lambda, XX_eigen, Lpart, theta_init_m, tol, max_iter, verbose))
  }
  
  if (model == "logistic") {
    thetahat = as.numeric(rnc_lasso_logistic_fit(y_m, X, gamma_n*Ln, lambda, theta_init_m, tol, max_iter, verbose))
  }
  
  if (model == "linear") {
    thetahat = as.numeric(rnc_lasso_linear_fit(y_m, X, gamma_n*Ln, lambda, theta_init_m, tol, max_iter, verbose))
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

rnc.lasso.predict = function(funkout, Xtest, Ln_full, train_ind) {
  
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
rnc.lasso.cv.single = function(y, X, Ln, reg_grid, train_ind, model = c("linear", "logistic", "poisson")) {
  
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
      rnc.m = rnc.lasso.fit(ytrain, Xtrain, Lntrain, reg_i, model)
      theta_ws = c(rnc.m$alphahat, rnc.m$betahat)
    }
    else {
      rnc.m = rnc.lasso.fit(ytrain, Xtrain, Lntrain, reg_i, model, theta_init = theta_ws)
    }
    pred = rnc.lasso.predict(rnc.m, Xtest, Ln, train_ind)

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
rnc.lasso.cv.Kfold = function(y, X, Ln, reg_grid, foldind, model = c("linear", "logistic", "poisson")) {
  n = nrow(X)
  nfolds = max(foldind)
  cv.ests = matrix(nrow = nrow(reg_grid), ncol = nfolds)
  for (k in 1:nfolds) {
    train_ind = which(foldind != k)
    cv.ests[,k] = rnc.lasso.cv.single(y, X, Ln, reg_grid, train_ind, model)
  }
  return(rowMeans(cv.ests))
}

rnc.lasso.cv.coorddesc = function(y, X, Ln, param_list, model = c("linear", "logistic", "poisson"), nfolds = 5, verb.cv = FALSE) {
  # coordinate descent with warm starts
  n = nrow(X)
  foldind = sample(nfolds, size = n, replace = TRUE) # change to be exactly n/K equal?
  
  gamma_n_grid = sort(param_list[[1]], decreasing = TRUE)
  lambda_grid = sort(param_list[[2]], decreasing = TRUE)
  
  gamma_n = gamma_n_grid[1]
  lambda = lambda_grid[1]
  reg_params_curr = c(gamma_n, lambda)
  cvg = FALSE
  iter = 0
  while (!cvg) {
    # Optimize gamma_p
    reg_grid = cbind(gamma_n_grid, lambda)
    cvKerr = rnc.lasso.cv.Kfold(y, X, Ln, reg_grid, foldind, model)
    gamma_n = gamma_n_grid[which.min(cvKerr)]
    if (verb.cv) print("# Optimized gamma_n #")
    # Optimize lambda
    reg_grid = cbind(gamma_n, lambda_grid)
    cvKerr = rnc.lasso.cv.Kfold(y, X, Ln, reg_grid, foldind, model)
    lambda = lambda_grid[which.min(cvKerr)]
    if (verb.cv) print("# Optimized lambda #")
    # Check convergence
    iter = iter + 1
    reg_params_new = c(gamma_n, lambda)
    step = sum(abs(reg_params_new-reg_params_curr))
    if (verb.cv) {
      print(paste0("### Finished coordinate descent iter ", iter, " ###"))
      print(paste0("Old params: ", paste0(reg_params_curr, collapse = " ")))
      print(paste0("New params: ", paste0(reg_params_new, collapse = " ")))
    }
    reg_params_curr = reg_params_new
    cvg = (step < 1e-04 || iter >= 10)
  }
  minerr = rnc.lasso.cv.Kfold(y, X, Ln, t(as.matrix(reg_params_curr)), foldind, model)
  return(list(reg_params_curr, minerr))  
}

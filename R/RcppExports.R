# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

SMOOTH <- function() {
    .Call(`_glmfunk_SMOOTH`)
}

loglike_bernoulli <- function(Y, P) {
    .Call(`_glmfunk_loglike_bernoulli`, Y, P)
}

logit_p <- function(eta) {
    .Call(`_glmfunk_logit_p`, eta)
}

proj_one <- function(a) {
    .Call(`_glmfunk_proj_one`, a)
}

soft_thresh <- function(a, q) {
    .Call(`_glmfunk_soft_thresh`, a, q)
}

prox_step <- function(theta, qx, n, p) {
    .Call(`_glmfunk_prox_step`, theta, qx, n, p)
}

proj_alpha <- function(Lp, w, mu, k) {
    .Call(`_glmfunk_proj_alpha`, Lp, w, mu, k)
}

prox_step_grace <- function(beta_tmp, qx, p) {
    .Call(`_glmfunk_prox_step_grace`, beta_tmp, qx, p)
}

logit_grace_l1_opt <- function(Y, X, Lp, lambda, mu, L, beta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_logit_grace_l1_opt`, Y, X, Lp, lambda, mu, L, beta_init, tol, max_iter, verbose)
}

grace_l1_logistic_fit <- function(Y, X, Lp, lambda, mu, L, beta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_grace_l1_logistic_fit`, Y, X, Lp, lambda, mu, L, beta_init, tol, max_iter, verbose)
}

linear_grace_l1_opt <- function(Y, X, Lp, lambda, mu, L, beta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_linear_grace_l1_opt`, Y, X, Lp, lambda, mu, L, beta_init, tol, max_iter, verbose)
}

grace_l1_linear_fit <- function(Y, X, Lp, lambda, mu, L, beta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_grace_l1_linear_fit`, Y, X, Lp, lambda, mu, L, beta_init, tol, max_iter, verbose)
}

logit_grace_l2_opt <- function(Y, X, Lp, lambda, beta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_logit_grace_l2_opt`, Y, X, Lp, lambda, beta_init, tol, max_iter, verbose)
}

grace_l2_logistic_fit <- function(Y, X, Lp, lambda, beta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_grace_l2_logistic_fit`, Y, X, Lp, lambda, beta_init, tol, max_iter, verbose)
}

linear_grace_l2_opt <- function(Y, X, Lp, lambda, beta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_linear_grace_l2_opt`, Y, X, Lp, lambda, beta_init, tol, max_iter, verbose)
}

grace_l2_linear_fit <- function(Y, X, Lp, lambda, beta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_grace_l2_linear_fit`, Y, X, Lp, lambda, beta_init, tol, max_iter, verbose)
}

poisson_funk_l1_opt <- function(Y, X, Ln, Lp, lambda, mu, XX_eigen, Lpart, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_poisson_funk_l1_opt`, Y, X, Ln, Lp, lambda, mu, XX_eigen, Lpart, theta_init, tol, max_iter, verbose)
}

funk_l1_poisson_fit <- function(Y, X, Ln, Lp, lambda, mu, XX_eigen, Lpart, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_funk_l1_poisson_fit`, Y, X, Ln, Lp, lambda, mu, XX_eigen, Lpart, theta_init, tol, max_iter, verbose)
}

logit_funk_l1_opt <- function(Y, X, Ln, Lp, lambda, mu, L, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_logit_funk_l1_opt`, Y, X, Ln, Lp, lambda, mu, L, theta_init, tol, max_iter, verbose)
}

funk_l1_logistic_fit <- function(Y, X, Ln, Lp, lambda, mu, L, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_funk_l1_logistic_fit`, Y, X, Ln, Lp, lambda, mu, L, theta_init, tol, max_iter, verbose)
}

linear_funk_l1_opt <- function(Y, X, Ln, Lp, lambda, mu, L, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_linear_funk_l1_opt`, Y, X, Ln, Lp, lambda, mu, L, theta_init, tol, max_iter, verbose)
}

funk_l1_linear_fit <- function(Y, X, Ln, Lp, lambda, mu, L, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_funk_l1_linear_fit`, Y, X, Ln, Lp, lambda, mu, L, theta_init, tol, max_iter, verbose)
}

poisson_funk_l2_opt <- function(Y, X, L, lambda, XX_eigen, Lpart, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_poisson_funk_l2_opt`, Y, X, L, lambda, XX_eigen, Lpart, theta_init, tol, max_iter, verbose)
}

funk_l2_poisson_fit <- function(Y, X, Ln, Lp, lambda, XX_eigen, Lpart, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_funk_l2_poisson_fit`, Y, X, Ln, Lp, lambda, XX_eigen, Lpart, theta_init, tol, max_iter, verbose)
}

logit_funk_l2_opt <- function(Y, X, L, lambda, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_logit_funk_l2_opt`, Y, X, L, lambda, theta_init, tol, max_iter, verbose)
}

funk_l2_logistic_fit <- function(Y, X, Ln, Lp, lambda, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_funk_l2_logistic_fit`, Y, X, Ln, Lp, lambda, theta_init, tol, max_iter, verbose)
}

linear_funk_l2_opt <- function(Y, X, L, lambda, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_linear_funk_l2_opt`, Y, X, L, lambda, theta_init, tol, max_iter, verbose)
}

funk_l2_linear_fit <- function(Y, X, Ln, Lp, lambda, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_funk_l2_linear_fit`, Y, X, Ln, Lp, lambda, theta_init, tol, max_iter, verbose)
}

poisson_rnc_lasso_opt <- function(Y, X, L, lambda, XX_eigen, Lpart, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_poisson_rnc_lasso_opt`, Y, X, L, lambda, XX_eigen, Lpart, theta_init, tol, max_iter, verbose)
}

rnc_lasso_poisson_fit <- function(Y, X, Ln, lambda, XX_eigen, Lpart, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_rnc_lasso_poisson_fit`, Y, X, Ln, lambda, XX_eigen, Lpart, theta_init, tol, max_iter, verbose)
}

logit_rnc_lasso_opt <- function(Y, X, L, lambda, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_logit_rnc_lasso_opt`, Y, X, L, lambda, theta_init, tol, max_iter, verbose)
}

rnc_lasso_logistic_fit <- function(Y, X, Ln, lambda, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_rnc_lasso_logistic_fit`, Y, X, Ln, lambda, theta_init, tol, max_iter, verbose)
}

linear_rnc_lasso_opt <- function(Y, X, L, lambda, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_linear_rnc_lasso_opt`, Y, X, L, lambda, theta_init, tol, max_iter, verbose)
}

rnc_lasso_linear_fit <- function(Y, X, Ln, lambda, theta_init, tol, max_iter, verbose) {
    .Call(`_glmfunk_rnc_lasso_linear_fit`, Y, X, Ln, lambda, theta_init, tol, max_iter, verbose)
}


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double SMOOTH() {
  return 0.00000001;
}

// logistic loss given y and fitted p
// [[Rcpp::export]]
double loglike_bernoulli(arma::mat Y, arma::mat P){
    double result = 0;
    int n = Y.n_rows;
    for (int i = 0; i<n; i++) {
        if (P(i,0) < SMOOTH()) {
            P(i,0) = SMOOTH();
        }
        if (P(i,0)>1-SMOOTH()) {
            P(i,0) = 1-SMOOTH();
        }
        result += Y(i,0)*log(P(i,0)) + (1-Y(i,0))*log(1-P(i,0));
    }
    return result;
}

// fitted p given \eta := X\beta
// [[Rcpp::export]]
arma::mat logit_p(arma::mat eta){
    int n = eta.n_rows;
    for (int i = 0; i<n; i++) {
        if (eta(i,0) > 12) {
            eta(i,0) = 12;
        }
        else if (eta(i,0) < -12) {
            eta(i,0) = -12;
        }
    }
    arma::mat expeta = exp(eta);
    arma::mat ones = arma::ones(n,1);
    arma::mat result = expeta/(expeta+ones);
    return result;
}

// project to l_infty ball with radius 1
// [[Rcpp::export]]
double proj_one(double a){
  if (fabs(a) < 1) {
    return a;
  }
  else if (a > 1) {
    return 1;
  }
  else {
    return (-1);
  }
}

// soft threshold a number at q
// [[Rcpp::export]]
double soft_thresh(double a, double q){
  if (a == 0) {
    return 0;
  }
  else if (a > 0) {
    return fmax(a - q, 0);
  }
  else {
    return (-1)*fmax(fabs(a) - q, 0);
  }
}

// prox theta: do nothing to alphas, soft-thresh betas at qx
// [[Rcpp::export]]
arma::mat prox_step(arma::mat theta, double qx, int n, int p){
    int nplusp = n + p;
    arma::vec theta_vec (nplusp);
    for(int i = 0; i < nplusp; ++i) {
      if (i < n) {
        theta_vec(i) = theta(i,0);
      }
      else {
        theta_vec(i) = soft_thresh(theta(i,0), qx);
      }
    }
    arma::mat theta_out=arma::mat(theta_vec);
    return theta_out;
}

// project everything
// [[Rcpp::export]]
arma::mat proj_alpha(arma::mat Lp, arma::mat w, double mu, int k){
    arma::mat s = (1/mu)*(Lp * w);    
    arma::vec w_vec (k);
    for(int i = 0; i < k; ++i) {
        w_vec(i) = proj_one(s(i,0));
    }
    arma::mat w_out = arma::mat(w_vec);
    return w_out;
}

// soft thresh all elements in beta_tmp at qx
// [[Rcpp::export]]
arma::mat prox_step_grace(arma::mat beta_tmp, double qx, int p){
    arma::vec beta_vec (p);
    for(int i = 0; i < p; ++i) {
        beta_vec(i) = soft_thresh(beta_tmp(i,0), qx);
    }
    arma::mat beta_out = arma::mat(beta_vec);
    return beta_out;
}

// function to solve logistic regression with l1 grace penalty
// [[Rcpp::export]]
arma::mat logit_grace_l1_opt(arma::mat Y, arma::mat X, arma::mat Lp, double lambda, double mu, double L, arma::mat beta_init, double tol, int max_iter, bool verbose){

    float n = X.n_rows;
    int p = X.n_cols;
    int k = Lp.n_rows;
    double st = lambda / L;
    
    int iter = 0;
    double err = 0;
    arma::mat beta_curr = beta_init;
    arma::mat w = beta_curr;
    double t_curr = 1;
    double t_prev = 1;
    bool converge = false;
    while (!converge) {
        arma::mat eta = X * w;
        arma::mat P = logit_p(eta);
        arma::mat residual = (1/n)*(P - Y);
        
        arma::mat alpha = proj_alpha(Lp, w, mu, k);
        arma::mat gradient1 = X.t()*residual;
        arma::mat gradient2 = Lp.t()*alpha;
        arma::mat gradient = gradient1 + gradient2;
        
        arma::mat beta_tmp = w - (1/L)*gradient;
        arma::mat beta_new = prox_step_grace(beta_tmp, st, p);
        
        t_prev = t_curr;
        t_curr = 2.0/(iter + 3.0);
        arma::mat wnew = beta_new + (1/t_prev - 1)*t_curr*(beta_new - beta_curr);

        err = arma::abs(beta_new - beta_curr).max();
        beta_curr = beta_new;
        w = wnew;

        if (err < tol & iter > 2000) {
            converge = true;
        }
        if(iter == max_iter){
            if (verbose) {
                Rcout << "Maximum iteration reached before convergence!" << std::endl;
            }

            break;
        }
        if (verbose) {
            Rcout << "Finished iteration " << iter << " with relative gap " << err << std::endl;
        }
        iter += 1;
    }
    return beta_curr;
}


// penalty parameters should be multiplied into both laplacians; lambda corresponds to lasso penalty
// [[Rcpp::export]]
extern "C" SEXP grace_l1_logistic_fit(SEXP Y, SEXP X, SEXP Lp, SEXP lambda, SEXP mu, SEXP L, SEXP beta_init, SEXP tol, SEXP max_iter, SEXP verbose){
    arma::mat Xmat=as<arma::mat>(X);
    arma::mat Lpdsmat=as<arma::mat>(Lp);
    arma::mat Ymat=as<arma::mat>(Y);
    double lambda_num= as<double>(lambda);
    double L_num= as<double>(L);
    double mu_num= as<double>(mu);
    double tol_num= as<double>(tol);
    bool verbose_ind= as<bool>(verbose);
    int iter_max_num= as<int>(max_iter);
    arma::mat betamat=as<arma::mat>(beta_init);
    arma::mat result = logit_grace_l1_opt(Ymat, Xmat, Lpdsmat, lambda_num, mu_num, L_num, betamat, tol_num, iter_max_num, verbose_ind);
    return(wrap(result));
}

// function to solve linear regression with l1 grace penalty
// [[Rcpp::export]]
arma::mat linear_grace_l1_opt(arma::mat Y, arma::mat X, arma::mat Lp, double lambda, double mu, double L, arma::mat beta_init, double tol, int max_iter, bool verbose){

    float n = X.n_rows;
    int p = X.n_cols;
    int k = Lp.n_rows;
    double st = lambda / L;
    
    int iter = 0;
    double err = 0;
    arma::mat beta_curr = beta_init;
    arma::mat w = beta_curr;
    double t_curr = 1;
    double t_prev = 1;
    bool converge = false;
    while (!converge) {
        arma::mat eta = X * w;
        arma::mat residual = (1/n)*(eta - Y);
        
        arma::mat alpha = proj_alpha(Lp, w, mu, k);
        arma::mat gradient1 = X.t()*residual;
        arma::mat gradient2 = Lp.t()*alpha;
        arma::mat gradient = gradient1 + gradient2;
        
        arma::mat beta_tmp = w - (1/L)*gradient;
        arma::mat beta_new = prox_step_grace(beta_tmp, st, p);
        
        t_prev = t_curr;
        t_curr = 2.0/(iter + 3.0);
        arma::mat wnew = beta_new + (1/t_prev - 1)*t_curr*(beta_new - beta_curr);

        err = arma::abs(beta_new - beta_curr).max();
        beta_curr = beta_new;
        w = wnew;

        if (err < tol & iter > 2000) {
            converge = true;
        }
        if(iter == max_iter){
            if (verbose) {
                Rcout << "Maximum iteration reached before convergence!" << std::endl;
            }

            break;
        }
        if (verbose) {
            Rcout << "Finished iteration " << iter << " with relative gap " << err << std::endl;
        }
        iter += 1;
    }
    return beta_curr;
}


// penalty parameters should be multiplied into both laplacians; lambda corresponds to lasso penalty
// [[Rcpp::export]]
extern "C" SEXP grace_l1_linear_fit(SEXP Y, SEXP X, SEXP Lp, SEXP lambda, SEXP mu, SEXP L, SEXP beta_init, SEXP tol, SEXP max_iter, SEXP verbose){
    arma::mat Xmat=as<arma::mat>(X);
    arma::mat Lpdsmat=as<arma::mat>(Lp);
    arma::mat Ymat=as<arma::mat>(Y);
    double lambda_num= as<double>(lambda);
    double L_num= as<double>(L);
    double mu_num= as<double>(mu);
    double tol_num= as<double>(tol);
    bool verbose_ind= as<bool>(verbose);
    int iter_max_num= as<int>(max_iter);
    arma::mat betamat=as<arma::mat>(beta_init);
    arma::mat result = linear_grace_l1_opt(Ymat, Xmat, Lpdsmat, lambda_num, mu_num, L_num, betamat, tol_num, iter_max_num, verbose_ind);
    return(wrap(result));
}


// function to solve logistic regression with l2 grace penalty
// [[Rcpp::export]]
arma::mat logit_grace_l2_opt(arma::mat Y, arma::mat X, arma::sp_mat Lp, double lambda, arma::mat beta_init, double tol, int max_iter, bool verbose){

    int p = X.n_cols;
    int iter = 0;
    double err = 0;
    arma::mat beta_old = beta_init;
    
    bool converge = false;
    while (!converge) {
        iter += 1;
        arma::mat eta = X * beta_old;
        arma::mat P = logit_p(eta);
        arma::mat residual = P - Y;
        arma::mat gradient = X.t()*residual + Lp*beta_old;
        arma::mat beta_tmp = beta_old - 0.001*gradient;
        arma::mat beta_new = prox_step_grace(beta_tmp, lambda, p);
        err = arma::norm(beta_new - beta_old,2)/(arma::norm(beta_old,2)+SMOOTH());
        if (err < tol) {
            converge = true;
        }
        beta_old = beta_new;
        if(iter == max_iter){
            if (verbose) {
                Rcout << "Maximum iteration reached before convergence!" << std::endl;
            }

            break;
        }
        if (verbose) {
            Rcout << "Finished iteration " << iter << " with relative gap " << err << std::endl;
        }
    }
    return beta_old;
}


// penalty parameters should be multiplied into both laplacians; lambda corresponds to lasso penalty
// [[Rcpp::export]]
extern "C" SEXP grace_l2_logistic_fit(SEXP Y, SEXP X, SEXP Lp, SEXP lambda, SEXP beta_init, SEXP tol, SEXP max_iter, SEXP verbose){
    arma::mat Xmat=as<arma::mat>(X);
    arma::mat Lpdsmat=as<arma::mat>(Lp);
    arma::sp_mat Lp_spmat=arma::sp_mat(Lpdsmat);
    arma::mat Ymat=as<arma::mat>(Y);
    double lambda_num= as<double>(lambda);
    double tol_num= as<double>(tol);
    bool verbose_ind= as<bool>(verbose);
    int iter_max_num= as<int>(max_iter);
    arma::mat betamat=as<arma::mat>(beta_init);
    arma::mat result = logit_grace_l2_opt(Ymat, Xmat, Lp_spmat, lambda_num, betamat, tol_num, iter_max_num, verbose_ind);
    return(wrap(result));
}


// function to solve logistic regression with l2 grace penalty
// [[Rcpp::export]]
arma::mat linear_grace_l2_opt(arma::mat Y, arma::mat X, arma::sp_mat Lp, double lambda, arma::mat beta_init, double tol, int max_iter, bool verbose){

    int p = X.n_cols;
    int iter = 0;
    double err = 0;
    arma::mat beta_old = beta_init;
    
    bool converge = false;
    while (!converge) {
        iter += 1;
        arma::mat eta = X * beta_old;
        arma::mat residual = eta - Y;
        arma::mat gradient = X.t()*residual + Lp*beta_old;
        arma::mat beta_tmp = beta_old - 0.001*gradient;
        arma::mat beta_new = prox_step_grace(beta_tmp, lambda, p);
        err = arma::norm(beta_new - beta_old,2)/(arma::norm(beta_old,2)+SMOOTH());
        if (err < tol) {
            converge = true;
        }
        beta_old = beta_new;
        if(iter == max_iter){
            if (verbose) {
                Rcout << "Maximum iteration reached before convergence!" << std::endl;
            }

            break;
        }
        if (verbose) {
            Rcout << "Finished iteration " << iter << " with relative gap " << err << std::endl;
        }
    }
    return beta_old;
}


// penalty parameters should be multiplied into both laplacians; lambda corresponds to lasso penalty
// [[Rcpp::export]]
extern "C" SEXP grace_l2_linear_fit(SEXP Y, SEXP X, SEXP Lp, SEXP lambda, SEXP beta_init, SEXP tol, SEXP max_iter, SEXP verbose){
    arma::mat Xmat=as<arma::mat>(X);
    arma::mat Lpdsmat=as<arma::mat>(Lp);
    arma::sp_mat Lp_spmat=arma::sp_mat(Lpdsmat);
    arma::mat Ymat=as<arma::mat>(Y);
    double lambda_num= as<double>(lambda);
    double tol_num= as<double>(tol);
    bool verbose_ind= as<bool>(verbose);
    int iter_max_num= as<int>(max_iter);
    arma::mat betamat=as<arma::mat>(beta_init);
    arma::mat result = linear_grace_l2_opt(Ymat, Xmat, Lp_spmat, lambda_num, betamat, tol_num, iter_max_num, verbose_ind);
    return(wrap(result));
}


// function to solve glm-funk Poisson regression with l1 feature network smoothing
// [[Rcpp::export]]
arma::mat poisson_funk_l1_opt(arma::mat Y, arma::mat X, arma::mat Ln, arma::mat Lp, double lambda, double mu, double XX_eigen, double Lpart, arma::mat theta_init, double tol, int max_iter, bool verbose){

    float n = X.n_rows;
    int p = X.n_cols;
    int k = Lp.n_rows;
    int nint = X.n_rows;
    arma::mat I = arma::eye<arma::mat>(nint,nint);
    arma::mat X_tilde = arma::join_rows(I,X);
    
    int iter = 0;
    double err = 0;
    arma::mat theta_curr = theta_init;
    arma::mat w = theta_curr;
    double t_curr = 1;
    double t_prev = 1;
    bool converge = false;

    while (!converge) {
        arma::mat eta = X_tilde * w;
        arma::mat P = exp(eta);
        double L = Lpart + P.max() * XX_eigen;
        double st = lambda / L;
        arma::mat residual = (1/n)*(P - Y);        
        arma::mat gradient1 = X_tilde.t()*residual;

        arma::mat w_p_submat = w.submat(n, 0, n + p - 1, 0);
        arma::mat alpha = proj_alpha(Lp, w_p_submat, mu, k);        
        arma::mat w_n_submat = w.submat(0, 0, n - 1, 0);
        arma::mat gradient2 = arma::join_cols(Ln*w_n_submat, Lp.t()*alpha);
        arma::mat gradient = gradient1 + gradient2;
        
        arma::mat theta_tmp = w - (1/L)*gradient;
        arma::mat theta_new = prox_step(theta_tmp, st, nint, p);
        t_prev = t_curr;
        t_curr = 2.0/(iter + 3.0);
        arma::mat wnew = theta_new + (1/t_prev - 1)*t_curr*(theta_new - theta_curr);

        err = arma::abs(theta_new - theta_curr).max();
        theta_curr = theta_new;
        w = wnew;

        if (err < tol & iter > 2000) {
            converge = true;
        }
        if(iter == max_iter){
            if (verbose) {
                Rcout << "Maximum iteration reached before convergence!" << std::endl;
            }

            break;
        }
        if (verbose) {
            Rcout << "Finished iteration " << iter << " with relative gap " << err << std::endl;
        }
        iter += 1;
    }
    return theta_curr;
}

// wrapper function; penalty parameters should be multiplied into both laplacians; lambda corresponds to lasso penalty
// [[Rcpp::export]]
extern "C" SEXP funk_l1_poisson_fit(SEXP Y, SEXP X, SEXP Ln, SEXP Lp, SEXP lambda, SEXP mu, SEXP XX_eigen, SEXP Lpart, SEXP theta_init, SEXP tol, SEXP max_iter, SEXP verbose){
    arma::mat Xmat=as<arma::mat>(X);
    arma::mat Lndsmat=as<arma::mat>(Ln);
    arma::mat Lpdsmat=as<arma::mat>(Lp);
    arma::mat Ymat=as<arma::mat>(Y);
    double lambda_num= as<double>(lambda);
    double XX_eigen_num= as<double>(XX_eigen);
    double Lpart_num= as<double>(Lpart);
    double mu_num= as<double>(mu);
    double tol_num= as<double>(tol);
    bool verbose_ind= as<bool>(verbose);
    int iter_max_num= as<int>(max_iter);
    arma::mat thetamat=as<arma::mat>(theta_init);
    arma::mat result = poisson_funk_l1_opt(Ymat, Xmat, Lndsmat, Lpdsmat, lambda_num, mu_num, XX_eigen_num, Lpart_num, thetamat, tol_num, iter_max_num, verbose_ind);
    return(wrap(result));
}


// function to solve glm-funk logistic regression with l1 feature network smoothing
// [[Rcpp::export]]
arma::mat logit_funk_l1_opt(arma::mat Y, arma::mat X, arma::mat Ln, arma::mat Lp, double lambda, double mu, double L, arma::mat theta_init, double tol, int max_iter, bool verbose){

    float n = X.n_rows;
    int p = X.n_cols;
    int k = Lp.n_rows;
    double st = lambda / L;
    int nint = X.n_rows;
    arma::mat I = arma::eye<arma::mat>(nint,nint);
    arma::mat X_tilde = arma::join_rows(I,X);
    
    int iter = 0;
    double err = 0;
    arma::mat theta_curr = theta_init;
    arma::mat w = theta_curr;
    double t_curr = 1;
    double t_prev = 1;
    bool converge = false;

    while (!converge) {
        arma::mat eta = X_tilde * w;
        arma::mat P = logit_p(eta);
        arma::mat residual = (1/n)*(P - Y);        
        arma::mat gradient1 = X_tilde.t()*residual;

        arma::mat w_p_submat = w.submat(n, 0, n + p - 1, 0);
        arma::mat alpha = proj_alpha(Lp, w_p_submat, mu, k);        
        arma::mat w_n_submat = w.submat(0, 0, n - 1, 0);
        arma::mat gradient2 = arma::join_cols(Ln*w_n_submat, Lp.t()*alpha);
        arma::mat gradient = gradient1 + gradient2;
        
        arma::mat theta_tmp = w - (1/L)*gradient;
        arma::mat theta_new = prox_step(theta_tmp, st, nint, p);
        t_prev = t_curr;
        t_curr = 2.0/(iter + 3.0);
        arma::mat wnew = theta_new + (1/t_prev - 1)*t_curr*(theta_new - theta_curr);

        err = arma::abs(theta_new - theta_curr).max();
        theta_curr = theta_new;
        w = wnew;

        if (err < tol & iter > 2000) {
            converge = true;
        }
        if(iter == max_iter){
            if (verbose) {
                Rcout << "Maximum iteration reached before convergence!" << std::endl;
            }

            break;
        }
        if (verbose) {
            Rcout << "Finished iteration " << iter << " with relative gap " << err << std::endl;
        }
        iter += 1;
    }
    return theta_curr;
}

// wrapper function; penalty parameters should be multiplied into both laplacians; lambda corresponds to lasso penalty
// [[Rcpp::export]]
extern "C" SEXP funk_l1_logistic_fit(SEXP Y, SEXP X, SEXP Ln, SEXP Lp, SEXP lambda, SEXP mu, SEXP L, SEXP theta_init, SEXP tol, SEXP max_iter, SEXP verbose){
    arma::mat Xmat=as<arma::mat>(X);
    arma::mat Lndsmat=as<arma::mat>(Ln);
    arma::mat Lpdsmat=as<arma::mat>(Lp);
    arma::mat Ymat=as<arma::mat>(Y);
    double lambda_num= as<double>(lambda);
    double L_num= as<double>(L);
    double mu_num= as<double>(mu);
    double tol_num= as<double>(tol);
    bool verbose_ind= as<bool>(verbose);
    int iter_max_num= as<int>(max_iter);
    arma::mat thetamat=as<arma::mat>(theta_init);
    arma::mat result = logit_funk_l1_opt(Ymat, Xmat, Lndsmat, Lpdsmat, lambda_num, mu_num, L_num, thetamat, tol_num, iter_max_num, verbose_ind);
    return(wrap(result));
}


// function to solve glm-funk linear regression with l1 feature network smoothing
// [[Rcpp::export]]
arma::mat linear_funk_l1_opt(arma::mat Y, arma::mat X, arma::mat Ln, arma::mat Lp, double lambda, double mu, double L, arma::mat theta_init, double tol, int max_iter, bool verbose){

    float n = X.n_rows;
    int p = X.n_cols;
    int k = Lp.n_rows;
    double st = lambda / L;
    int nint = X.n_rows;
    arma::mat I = arma::eye<arma::mat>(nint,nint);
    arma::mat X_tilde = arma::join_rows(I,X);
    
    int iter = 0;
    double err = 0;
    arma::mat theta_curr = theta_init;
    arma::mat w = theta_curr;
    double t_curr = 1;
    double t_prev = 1;
    bool converge = false;

    while (!converge) {
        arma::mat eta = X_tilde * w;
        arma::mat residual = (1/n)*(eta - Y);        
        arma::mat gradient1 = X_tilde.t()*residual;

        arma::mat w_p_submat = w.submat(n, 0, n + p - 1, 0);
        arma::mat alpha = proj_alpha(Lp, w_p_submat, mu, k);        
        arma::mat w_n_submat = w.submat(0, 0, n - 1, 0);
        arma::mat gradient2 = arma::join_cols(Ln*w_n_submat, Lp.t()*alpha);
        arma::mat gradient = gradient1 + gradient2;
        
        arma::mat theta_tmp = w - (1/L)*gradient;
        arma::mat theta_new = prox_step(theta_tmp, st, nint, p);
        t_prev = t_curr;
        t_curr = 2.0/(iter + 3.0);
        arma::mat wnew = theta_new + (1/t_prev - 1)*t_curr*(theta_new - theta_curr);

        err = arma::abs(theta_new - theta_curr).max();
        theta_curr = theta_new;
        w = wnew;

        if (err < tol & iter > 2000) {
            converge = true;
        }
        if(iter == max_iter){
            if (verbose) {
                Rcout << "Maximum iteration reached before convergence!" << std::endl;
            }

            break;
        }
        if (verbose) {
            Rcout << "Finished iteration " << iter << " with relative gap " << err << std::endl;
        }
        iter += 1;
    }
    return theta_curr;
}

// wrapper function; penalty parameters should be multiplied into both laplacians; lambda corresponds to lasso penalty
// [[Rcpp::export]]
extern "C" SEXP funk_l1_linear_fit(SEXP Y, SEXP X, SEXP Ln, SEXP Lp, SEXP lambda, SEXP mu, SEXP L, SEXP theta_init, SEXP tol, SEXP max_iter, SEXP verbose){
    arma::mat Xmat=as<arma::mat>(X);
    arma::mat Lndsmat=as<arma::mat>(Ln);
    arma::mat Lpdsmat=as<arma::mat>(Lp);
    arma::mat Ymat=as<arma::mat>(Y);
    double lambda_num= as<double>(lambda);
    double L_num= as<double>(L);
    double mu_num= as<double>(mu);
    double tol_num= as<double>(tol);
    bool verbose_ind= as<bool>(verbose);
    int iter_max_num= as<int>(max_iter);
    arma::mat thetamat=as<arma::mat>(theta_init);
    arma::mat result = linear_funk_l1_opt(Ymat, Xmat, Lndsmat, Lpdsmat, lambda_num, mu_num, L_num, thetamat, tol_num, iter_max_num, verbose_ind);
    return(wrap(result));
}


// function to solve glm-funk Poisson regression with l2 feature network smoothing
// [[Rcpp::export]]
arma::mat poisson_funk_l2_opt(arma::mat Y, arma::mat X, arma::sp_mat L, double lambda, double XX_eigen, double Lpart, arma::mat theta_init, double tol, int max_iter, bool verbose){

    int n = X.n_rows;
    int p = X.n_cols;
    int iter = 0;
    double err = 0;
    arma::mat I = arma::eye<arma::mat>(n,n);
    arma::mat X_tilde = arma::join_rows(I,X);
    arma::mat theta_old = theta_init;
    
    bool converge = false;
    while (!converge) {
        iter += 1;
        arma::mat eta = X_tilde * theta_old;        
        arma::mat P = exp(eta);
        double C_L = Lpart + P.max() * XX_eigen;
        arma::mat residual = P - Y;
        arma::mat gradient = X_tilde.t()*residual + L*theta_old;
        arma::mat theta_tmp = theta_old - (1 / C_L)*gradient;
        arma::mat theta_new = prox_step(theta_tmp, lambda, n, p);
        err = arma::norm(theta_new - theta_old,2)/(arma::norm(theta_old,2)+SMOOTH());
        if (err < tol) {
            converge = true;
        }
        theta_old = theta_new;
        if(iter == max_iter){
            if (verbose) {
                Rcout << "Maximum iteration reached before convergence!" << std::endl;
            }

            break;
        }
        if (verbose) {
            Rcout << "Finished iteration " << iter << " with relative gap " << err << std::endl;
        }
    }
    return theta_old;
}


// wrapper function; penalty parameters should be multiplied into both laplacians; lambda corresponds to lasso penalty
// [[Rcpp::export]]
extern "C" SEXP funk_l2_poisson_fit(SEXP Y, SEXP X, SEXP Ln, SEXP Lp, SEXP lambda, SEXP XX_eigen, SEXP Lpart, SEXP theta_init, SEXP tol, SEXP max_iter, SEXP verbose){
    arma::mat Xmat=as<arma::mat>(X);
    arma::mat Lndsmat=as<arma::mat>(Ln);
    arma::mat Lpdsmat=as<arma::mat>(Lp);
    int n = Xmat.n_rows;
    int p = Xmat.n_cols;
    arma::mat Omega = arma::zeros<arma::mat>(n+p,n+p);
    Omega(arma::span(0,n-1),arma::span(0,n-1)) = Lndsmat;
    Omega(arma::span(n,n+p-1),arma::span(n,n+p-1)) = Lpdsmat;
    arma::sp_mat Lspmat=arma::sp_mat(Omega);
    arma::mat Ymat=as<arma::mat>(Y);
    double lambda_num= as<double>(lambda);
    double XX_eigen_num= as<double>(XX_eigen);
    double Lpart_num= as<double>(Lpart);    
    double tol_num= as<double>(tol);
    bool verbose_ind= as<bool>(verbose);
    int iter_max_num= as<int>(max_iter);
    arma::mat thetamat=as<arma::mat>(theta_init);
    arma::mat result = poisson_funk_l2_opt(Ymat, Xmat, Lspmat, lambda_num, XX_eigen_num, Lpart_num, thetamat, tol_num, iter_max_num, verbose_ind);
    return(wrap(result));
}


// function to solve glm-funk logistic regression with l2 feature network smoothing
// [[Rcpp::export]]
arma::mat logit_funk_l2_opt(arma::mat Y, arma::mat X, arma::sp_mat L, double lambda, arma::mat theta_init, double tol, int max_iter, bool verbose){

    int n = X.n_rows;
    int p = X.n_cols;
    int iter = 0;
    double err = 0;
    arma::mat I = arma::eye<arma::mat>(n,n);
    arma::mat X_tilde = arma::join_rows(I,X);
    arma::mat theta_old = theta_init;
    
    bool converge = false;
    while (!converge) {
        iter += 1;
        arma::mat eta = X_tilde * theta_old;
        arma::mat P = logit_p(eta);
        arma::mat residual = P - Y;
        arma::mat gradient = X_tilde.t()*residual + L*theta_old;
        arma::mat theta_tmp = theta_old - 0.001*gradient;
        arma::mat theta_new = prox_step(theta_tmp, lambda, n, p);
        err = arma::norm(theta_new - theta_old,2)/(arma::norm(theta_old,2)+SMOOTH());
        if (err < tol) {
            converge = true;
        }
        theta_old = theta_new;
        if(iter == max_iter){
            if (verbose) {
                Rcout << "Maximum iteration reached before convergence!" << std::endl;
            }

            break;
        }
        if (verbose) {
            Rcout << "Finished iteration " << iter << " with relative gap " << err << std::endl;
        }
    }
    return theta_old;
}


// wrapper function; penalty parameters should be multiplied into both laplacians; lambda corresponds to lasso penalty
// [[Rcpp::export]]
extern "C" SEXP funk_l2_logistic_fit(SEXP Y, SEXP X, SEXP Ln, SEXP Lp, SEXP lambda, SEXP theta_init, SEXP tol, SEXP max_iter, SEXP verbose){
    arma::mat Xmat=as<arma::mat>(X);
    arma::mat Lndsmat=as<arma::mat>(Ln);
    arma::mat Lpdsmat=as<arma::mat>(Lp);
    int n = Xmat.n_rows;
    int p = Xmat.n_cols;
    arma::mat Omega = arma::zeros<arma::mat>(n+p,n+p);
    Omega(arma::span(0,n-1),arma::span(0,n-1)) = Lndsmat;
    Omega(arma::span(n,n+p-1),arma::span(n,n+p-1)) = Lpdsmat;
    arma::sp_mat Lspmat=arma::sp_mat(Omega);
    arma::mat Ymat=as<arma::mat>(Y);
    double lambda_num= as<double>(lambda);
    double tol_num= as<double>(tol);
    bool verbose_ind= as<bool>(verbose);
    int iter_max_num= as<int>(max_iter);
    arma::mat thetamat=as<arma::mat>(theta_init);
    arma::mat result = logit_funk_l2_opt(Ymat, Xmat, Lspmat, lambda_num, thetamat, tol_num, iter_max_num, verbose_ind);
    return(wrap(result));
}


// function to solve glm-funk linear regression with l2 feature network smoothing
// [[Rcpp::export]]
arma::mat linear_funk_l2_opt(arma::mat Y, arma::mat X, arma::sp_mat L, double lambda, arma::mat theta_init, double tol, int max_iter, bool verbose){

    int n = X.n_rows;
    int p = X.n_cols;
    int iter = 0;
    double err = 0;
    arma::mat I = arma::eye<arma::mat>(n,n);
    arma::mat X_tilde = arma::join_rows(I,X);
    arma::mat theta_old = theta_init;
    
    bool converge = false;
    while (!converge) {
        iter += 1;
        arma::mat eta = X_tilde * theta_old;
        arma::mat residual = eta - Y;
        arma::mat gradient = X_tilde.t()*residual + L*theta_old;
        arma::mat theta_tmp = theta_old - 0.001*gradient;
        arma::mat theta_new = prox_step(theta_tmp, lambda, n, p);
        err = arma::norm(theta_new - theta_old,2)/(arma::norm(theta_old,2)+SMOOTH());
        if (err < tol) {
            converge = true;
        }
        theta_old = theta_new;
        if(iter == max_iter){
            if (verbose) {
                Rcout << "Maximum iteration reached before convergence!" << std::endl;
            }

            break;
        }
        if (verbose) {
            Rcout << "Finish iteration " << iter << " with relative gap " << err << std::endl;
        }
    }
    return theta_old;
}


// wrapper function; penalty parameters should be multiplied into both laplacians; lambda corresponds to lasso penalty
// [[Rcpp::export]]
extern "C" SEXP funk_l2_linear_fit(SEXP Y, SEXP X, SEXP Ln, SEXP Lp, SEXP lambda, SEXP theta_init, SEXP tol, SEXP max_iter, SEXP verbose){
    arma::mat Xmat=as<arma::mat>(X);
    arma::mat Lndsmat=as<arma::mat>(Ln);
    arma::mat Lpdsmat=as<arma::mat>(Lp);
    int n = Xmat.n_rows;
    int p = Xmat.n_cols;
    arma::mat Omega = arma::zeros<arma::mat>(n+p,n+p);
    Omega(arma::span(0,n-1),arma::span(0,n-1)) = Lndsmat;
    Omega(arma::span(n,n+p-1),arma::span(n,n+p-1)) = Lpdsmat;
    arma::sp_mat Lspmat=arma::sp_mat(Omega);
    arma::mat Ymat=as<arma::mat>(Y);
    double lambda_num= as<double>(lambda);
    double tol_num= as<double>(tol);
    bool verbose_ind= as<bool>(verbose);
    int iter_max_num= as<int>(max_iter);
    arma::mat thetamat=as<arma::mat>(theta_init);
    arma::mat result = linear_funk_l2_opt(Ymat, Xmat, Lspmat, lambda_num, thetamat, tol_num, iter_max_num, verbose_ind);
    return(wrap(result));
}


// function to solve Poisson regression with network cohesion, when covariates are provided
// [[Rcpp::export]]
arma::mat poisson_rnc_lasso_opt(arma::mat Y, arma::mat X, arma::mat L, double lambda, double XX_eigen, double Lpart, arma::mat theta_init, double tol, int max_iter, bool verbose){

    int n = X.n_rows;
    int p = X.n_cols;
    int iter = 0;
    double err = 0;
    arma::mat I = arma::eye<arma::mat>(n,n);
    arma::mat X_tilde = arma::join_rows(I,X);
    arma::mat theta_old = theta_init;
    
    bool converge = false;
    while (!converge) {      
        iter += 1;
        arma::mat eta = X_tilde * theta_old;
        arma::mat P = exp(eta);
        double C_L = Lpart + P.max() * XX_eigen;
        arma::mat residual = P - Y;
        arma::mat gradient = X_tilde.t()*residual + L*theta_old;
        arma::mat theta_tmp = theta_old - (1 / C_L)*gradient;
        arma::mat theta_new = prox_step(theta_tmp, lambda, n, p);
        err = arma::norm(theta_new - theta_old,2)/(arma::norm(theta_old,2)+SMOOTH());
        if (err < tol) {
            converge = true;
        }
        theta_old = theta_new;
        if(iter == max_iter){
            if (verbose) {
                Rcout << "Maximum iteration reached before convergence!" << std::endl;
            }

            break;
        }
        if (verbose) {
            Rcout << "Finished iteration " << iter << " with relative gap " << err << std::endl;
        }
    }
    return theta_old;
}

// penalty parameters should be multiplied into both laplacians; lambda corresponds to lasso penalty
// [[Rcpp::export]]
extern "C" SEXP rnc_lasso_poisson_fit(SEXP Y, SEXP X, SEXP Ln, SEXP lambda, SEXP XX_eigen, SEXP Lpart, SEXP theta_init, SEXP tol, SEXP max_iter, SEXP verbose){
    arma::mat Xmat=as<arma::mat>(X);
    arma::mat Lndsmat=as<arma::mat>(Ln);
    int n = Xmat.n_rows;
    int p = Xmat.n_cols;
    arma::mat Omega = arma::zeros<arma::mat>(n+p,n+p);
    Omega(arma::span(0,n-1),arma::span(0,n-1)) = Lndsmat;    
    arma::mat Ymat=as<arma::mat>(Y);
    double lambda_num= as<double>(lambda);
    double XX_eigen_num= as<double>(XX_eigen);
    double Lpart_num= as<double>(Lpart);    
    double tol_num= as<double>(tol);
    bool verbose_ind= as<bool>(verbose);
    int iter_max_num= as<int>(max_iter);
    arma::mat thetamat=as<arma::mat>(theta_init);
    arma::mat result = poisson_rnc_lasso_opt(Ymat, Xmat, Omega, lambda_num, XX_eigen_num, Lpart_num, thetamat, tol_num, iter_max_num, verbose_ind);
    return(wrap(result));
}



// function to solve logistic regression with network cohesion, when covariates are provided
// [[Rcpp::export]]
arma::mat logit_rnc_lasso_opt(arma::mat Y, arma::mat X, arma::mat L, double lambda, arma::mat theta_init, double tol, int max_iter, bool verbose){

    int n = X.n_rows;
    int p = X.n_cols;
    int iter = 0;
    double err = 0;
    arma::mat I = arma::eye<arma::mat>(n,n);
    arma::mat X_tilde = arma::join_rows(I,X);
    arma::mat theta_old = theta_init;
    
    bool converge = false;
    while (!converge) {
        iter += 1;
        arma::mat eta = X_tilde * theta_old;
        arma::mat P = logit_p(eta);
        arma::mat residual = P - Y;
        arma::mat gradient = X_tilde.t()*residual + L*theta_old;
        arma::mat theta_tmp = theta_old - 0.001*gradient;
        arma::mat theta_new = prox_step(theta_tmp, lambda, n, p);
        err = arma::norm(theta_new - theta_old,2)/(arma::norm(theta_old,2)+SMOOTH());
        if (err < tol) {
            converge = true;
        }
        theta_old = theta_new;
        if(iter == max_iter){
            if (verbose) {
                Rcout << "Maximum iteration reached before convergence!" << std::endl;
            }

            break;
        }
        if (verbose) {
            Rcout << "Finished iteration " << iter << " with relative gap " << err << std::endl;
        }
    }
    return theta_old;
}


// penalty parameters should be multiplied into both laplacians; lambda corresponds to lasso penalty
// [[Rcpp::export]]
extern "C" SEXP rnc_lasso_logistic_fit(SEXP Y, SEXP X, SEXP Ln, SEXP lambda, SEXP theta_init, SEXP tol, SEXP max_iter, SEXP verbose){
    arma::mat Xmat=as<arma::mat>(X);
    arma::mat Lndsmat=as<arma::mat>(Ln);
    int n = Xmat.n_rows;
    int p = Xmat.n_cols;
    arma::mat Omega = arma::zeros<arma::mat>(n+p,n+p);
    Omega(arma::span(0,n-1),arma::span(0,n-1)) = Lndsmat;    
    arma::mat Ymat=as<arma::mat>(Y);
    double lambda_num= as<double>(lambda);
    double tol_num= as<double>(tol);
    bool verbose_ind= as<bool>(verbose);
    int iter_max_num= as<int>(max_iter);
    arma::mat thetamat=as<arma::mat>(theta_init);
    arma::mat result = logit_rnc_lasso_opt(Ymat, Xmat, Omega, lambda_num, thetamat, tol_num, iter_max_num, verbose_ind);
    return(wrap(result));
}


// function to solve linear regression with network cohesion, when covariates are provided
// [[Rcpp::export]]
arma::mat linear_rnc_lasso_opt(arma::mat Y, arma::mat X, arma::mat L, double lambda, arma::mat theta_init, double tol, int max_iter, bool verbose){

    int n = X.n_rows;
    int p = X.n_cols;
    int iter = 0;
    double err = 0;
    arma::mat I = arma::eye<arma::mat>(n,n);
    arma::mat X_tilde = arma::join_rows(I,X);
    arma::mat theta_old = theta_init;
    
    bool converge = false;
    while (!converge) {
        iter += 1;
        arma::mat eta = X_tilde * theta_old;
        arma::mat residual = eta - Y;
        arma::mat gradient = X_tilde.t()*residual + L*theta_old;
        arma::mat theta_tmp = theta_old - 0.001*gradient;
        arma::mat theta_new = prox_step(theta_tmp, lambda, n, p);
        err = arma::norm(theta_new - theta_old,2)/(arma::norm(theta_old,2)+SMOOTH());
        if (err < tol) {
            converge = true;
        }
        theta_old = theta_new;
        if(iter == max_iter){
            if (verbose) {
                Rcout << "Maximum iteration reached before convergence!" << std::endl;
            }

            break;
        }
        if (verbose) {
            Rcout << "Finished iteration " << iter << " with relative gap " << err << std::endl;
        }
    }
    return theta_old;
}


// penalty parameters should be multiplied into both laplacians; lambda corresponds to lasso penalty
// [[Rcpp::export]]
extern "C" SEXP rnc_lasso_linear_fit(SEXP Y, SEXP X, SEXP Ln, SEXP lambda, SEXP theta_init, SEXP tol, SEXP max_iter, SEXP verbose){
    arma::mat Xmat=as<arma::mat>(X);
    arma::mat Lndsmat=as<arma::mat>(Ln);
    int n = Xmat.n_rows;
    int p = Xmat.n_cols;
    arma::mat Omega = arma::zeros<arma::mat>(n+p,n+p);
    Omega(arma::span(0,n-1),arma::span(0,n-1)) = Lndsmat;    
    arma::mat Ymat=as<arma::mat>(Y);
    double lambda_num= as<double>(lambda);
    double tol_num= as<double>(tol);
    bool verbose_ind= as<bool>(verbose);
    int iter_max_num= as<int>(max_iter);
    arma::mat thetamat=as<arma::mat>(theta_init);
    arma::mat result = linear_rnc_lasso_opt(Ymat, Xmat, Omega, lambda_num, thetamat, tol_num, iter_max_num, verbose_ind);
    return(wrap(result));
}


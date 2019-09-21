#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _glmfunk_funk_l1_linear_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_funk_l1_logistic_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_funk_l2_linear_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_funk_l2_logistic_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_grace_l1_linear_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_grace_l1_logistic_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_grace_l2_linear_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_grace_l2_logistic_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_linear_funk_l1_opt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_linear_funk_l2_opt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_linear_grace_l1_opt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_linear_grace_l2_opt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_linear_rnc_lasso_opt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_logit_funk_l1_opt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_logit_funk_l2_opt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_logit_grace_l1_opt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_logit_grace_l2_opt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_logit_p(SEXP);
extern SEXP _glmfunk_logit_rnc_lasso_opt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_loglike_bernoulli(SEXP, SEXP);
extern SEXP _glmfunk_proj_alpha(SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_proj_one(SEXP);
extern SEXP _glmfunk_prox_step(SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_prox_step_grace(SEXP, SEXP, SEXP);
extern SEXP _glmfunk_rnc_lasso_linear_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_rnc_lasso_logistic_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmfunk_SMOOTH();
extern SEXP _glmfunk_soft_thresh(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_glmfunk_funk_l1_linear_fit",     (DL_FUNC) &_glmfunk_funk_l1_linear_fit,     11},
    {"_glmfunk_funk_l1_logistic_fit",   (DL_FUNC) &_glmfunk_funk_l1_logistic_fit,   11},
    {"_glmfunk_funk_l2_linear_fit",     (DL_FUNC) &_glmfunk_funk_l2_linear_fit,      9},
    {"_glmfunk_funk_l2_logistic_fit",   (DL_FUNC) &_glmfunk_funk_l2_logistic_fit,    9},
    {"_glmfunk_grace_l1_linear_fit",    (DL_FUNC) &_glmfunk_grace_l1_linear_fit,    10},
    {"_glmfunk_grace_l1_logistic_fit",  (DL_FUNC) &_glmfunk_grace_l1_logistic_fit,  10},
    {"_glmfunk_grace_l2_linear_fit",    (DL_FUNC) &_glmfunk_grace_l2_linear_fit,     8},
    {"_glmfunk_grace_l2_logistic_fit",  (DL_FUNC) &_glmfunk_grace_l2_logistic_fit,   8},
    {"_glmfunk_linear_funk_l1_opt",     (DL_FUNC) &_glmfunk_linear_funk_l1_opt,     11},
    {"_glmfunk_linear_funk_l2_opt",     (DL_FUNC) &_glmfunk_linear_funk_l2_opt,      8},
    {"_glmfunk_linear_grace_l1_opt",    (DL_FUNC) &_glmfunk_linear_grace_l1_opt,    10},
    {"_glmfunk_linear_grace_l2_opt",    (DL_FUNC) &_glmfunk_linear_grace_l2_opt,     8},
    {"_glmfunk_linear_rnc_lasso_opt",   (DL_FUNC) &_glmfunk_linear_rnc_lasso_opt,    8},
    {"_glmfunk_logit_funk_l1_opt",      (DL_FUNC) &_glmfunk_logit_funk_l1_opt,      11},
    {"_glmfunk_logit_funk_l2_opt",      (DL_FUNC) &_glmfunk_logit_funk_l2_opt,       8},
    {"_glmfunk_logit_grace_l1_opt",     (DL_FUNC) &_glmfunk_logit_grace_l1_opt,     10},
    {"_glmfunk_logit_grace_l2_opt",     (DL_FUNC) &_glmfunk_logit_grace_l2_opt,      8},
    {"_glmfunk_logit_p",                (DL_FUNC) &_glmfunk_logit_p,                 1},
    {"_glmfunk_logit_rnc_lasso_opt",    (DL_FUNC) &_glmfunk_logit_rnc_lasso_opt,     8},
    {"_glmfunk_loglike_bernoulli",      (DL_FUNC) &_glmfunk_loglike_bernoulli,       2},
    {"_glmfunk_proj_alpha",             (DL_FUNC) &_glmfunk_proj_alpha,              4},
    {"_glmfunk_proj_one",               (DL_FUNC) &_glmfunk_proj_one,                1},
    {"_glmfunk_prox_step",              (DL_FUNC) &_glmfunk_prox_step,               4},
    {"_glmfunk_prox_step_grace",        (DL_FUNC) &_glmfunk_prox_step_grace,         3},
    {"_glmfunk_rnc_lasso_linear_fit",   (DL_FUNC) &_glmfunk_rnc_lasso_linear_fit,    8},
    {"_glmfunk_rnc_lasso_logistic_fit", (DL_FUNC) &_glmfunk_rnc_lasso_logistic_fit,  8},
    {"_glmfunk_SMOOTH",                 (DL_FUNC) &_glmfunk_SMOOTH,                  0},
    {"_glmfunk_soft_thresh",            (DL_FUNC) &_glmfunk_soft_thresh,             2},
    {NULL, NULL, 0}
};

void R_init_glmfunk(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

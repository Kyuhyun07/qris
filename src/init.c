#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _qrismb_Amat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _qrismb_isObj(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _qrismb_rev_isObj(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_qrismb_Amat",      (DL_FUNC) &_qrismb_Amat,      7},
    {"_qrismb_isObj",     (DL_FUNC) &_qrismb_isObj,     7},
    {"_qrismb_rev_isObj", (DL_FUNC) &_qrismb_rev_isObj, 8},
    {NULL, NULL, 0}
};

void R_init_qrismb(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

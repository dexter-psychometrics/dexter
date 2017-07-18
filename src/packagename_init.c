#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void E(void *, void *, void *, void *, void *, void *, void *, void *);
extern void ElSym(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void E_n(void *, void *, void *, void *, void *, void *, void *, void *);
extern void Escore(void *, void *, void *, void *, void *, void *, void *);
extern void H(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ittot_mat(void *, void *, void *, void *, void *, void *, void *, void *);
extern void ittot_mat0(void *, void *, void *, void *, void *, void *, void *, void *);
extern void meanElSym(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void PV(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sampleNRM(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"E",          (DL_FUNC) &E,           8},
  {"ElSym",      (DL_FUNC) &ElSym,       9},
  {"E_n",        (DL_FUNC) &E_n,         8},
  {"Escore",     (DL_FUNC) &Escore,      7},
  {"H",          (DL_FUNC) &H,           9},
  {"ittot_mat",  (DL_FUNC) &ittot_mat,   8},
  {"ittot_mat0", (DL_FUNC) &ittot_mat0,  8},
  {"meanElSym",  (DL_FUNC) &meanElSym,  10},
  {"PV",         (DL_FUNC) &PV,         13},
  {"sampleNRM",  (DL_FUNC) &sampleNRM,   8},
  {NULL, NULL, 0}
};

void R_init_dexter(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

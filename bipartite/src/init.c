#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* 
generated with
tools::package_native_routine_registration_skeleton(".")
*/


/* .C calls */
/* extern void bmn5(int*, int*, int*, double*, int*, int*, int*,  int*,  int*,  int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int*, int*); */
extern void identifyModules(int*, char**);

static const R_CMethodDef cMethods[ ] = {
/*    {"bmn5", (DL_FUNC) &bmn5, 21}, */
    {"identifyModules", (DL_FUNC) &identifyModules,  2 },
    {NULL, NULL, 0}
};

/* {INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP,  INTSXP,  INTSXP,  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP} */
/* {INTSXP, STRSXP}*/

void R_init_bipartite(DllInfo *dll)
{
    R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

void R_unload_bipartite(DllInfo *dll)
{
  /* resease resources */
}

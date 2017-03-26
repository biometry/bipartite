#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* 
generated with
tools::package_native_routine_registration_skeleton(".")
*/


/* .C calls */
extern void bmn5(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void identifyModules(void *, void *);

static const R_CMethodDef CallEntries[ ] = {
    {"bmn5",            (DL_FUNC) &bmn5,            21},
    {"identifyModules", (DL_FUNC) &identifyModules,  2},
    {NULL, NULL, 0}
};

void R_init_bipartite(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
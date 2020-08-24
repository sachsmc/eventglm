#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

    /* .C calls */
extern void loo_comprisk2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void loo_comprisk(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void loo_surv2(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void loo_surv(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"loo_comprisk2", (DL_FUNC) &loo_comprisk2, 11},
    {"loo_comprisk", (DL_FUNC) &loo_comprisk, 11},
    {"loo_surv2",     (DL_FUNC) &loo_surv2,      9},
    {"loo_surv",     (DL_FUNC) &loo_surv,      9},
    {NULL, NULL, 0}
};

void R_init_eventglm(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

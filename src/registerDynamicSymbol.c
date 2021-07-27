#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
  Check these declarations against the C/Fortran source code.
*/

  /* .C calls */
extern void cpChangePointDStat(void *, void *, void *, void *, void *);
extern void cpCopulaStats(void *, void *, void *, void *, void *, void *);
extern void cpCopulaStatsBucher(void *, void *, void *, void *, void *, void *);
extern void cpCopulaStatsMult(void *, void *, void *, void *, void *, void *);
extern void cpCopulaStatsMultBucherNonSeq(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void cpCopulaStatsMultBucherSeq(void *, void *, void *, void *, void *, void *, void *);
extern void empcdf(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"cpChangePointDStat",            (DL_FUNC) &cpChangePointDStat,            5},
  {"cpCopulaStats",                 (DL_FUNC) &cpCopulaStats,                 6},
  {"cpCopulaStatsBucher",           (DL_FUNC) &cpCopulaStatsBucher,           6},
  {"cpCopulaStatsMult",             (DL_FUNC) &cpCopulaStatsMult,             6},
  {"cpCopulaStatsMultBucherNonSeq", (DL_FUNC) &cpCopulaStatsMultBucherNonSeq, 9},
  {"cpCopulaStatsMultBucherSeq",    (DL_FUNC) &cpCopulaStatsMultBucherSeq,    7},
  {"empcdf",                        (DL_FUNC) &empcdf,                        6},
  {NULL, NULL, 0}
};

void R_init_changepointTests(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

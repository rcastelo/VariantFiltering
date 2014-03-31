#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include "VariantFiltering.h"

static R_CallMethodDef callMethods[] = {
  /* methods-WeightMatrix.c */
  {"scoss_read_wm", (DL_FUNC) &scoss_read_wm, 1},
  {"scoss_wm_score", (DL_FUNC) &scoss_wm_score, 3},
  {"scoss_wm_score_DNAStringSet", (DL_FUNC) &scoss_wm_score_DNAStringSet, 3},
  {"scoss_width_wm", (DL_FUNC) &scoss_width_wm, 1},
  {"scoss_conserved_positions_wm", (DL_FUNC) &scoss_conserved_positions_wm, 1},
  {"scoss_show_wm", (DL_FUNC) &scoss_show_wm, 1},
  {NULL, NULL, 0}
};

void
R_init_VariantFiltering(DllInfo* info) {

  R_registerRoutines(info, NULL, callMethods, NULL, NULL);

}

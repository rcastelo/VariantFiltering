#ifndef _VARIANTFILTERING_H_
#define _VARIANTFILTERING_H_

#include <Rdefines.h>

SEXP
scoss_read_wm(SEXP fnameR);

SEXP
scoss_wm_score(SEXP wmR, SEXP dnastringR, SEXP nscoR);

SEXP
scoss_wm_score_DNAStringSet(SEXP wmR, SEXP dnastringR, SEXP nscoR);

SEXP
scoss_width_wm(SEXP wmR);

SEXP
scoss_conserved_positions_wm(SEXP wmR);

void
scoss_show_wm(SEXP wmR);

#endif                                /* _VARIANTFILTERING_H_ */

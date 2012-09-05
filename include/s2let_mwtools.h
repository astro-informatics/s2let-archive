// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_MWTOOLS
#define S2LET_MWTOOLS

void s2let_mw_alm2map(complex double* f, const complex double* flm, int L);
void s2let_mw_map2alm(complex double* flm, const complex double* f, int L);
void s2let_mw_alm2map_real(double* f, const complex double* flm, int L);
void s2let_mw_map2alm_real(complex double* flm, const double* f, int L);

#endif

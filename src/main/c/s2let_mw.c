// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <complex.h>
#include <ssht.h>
#include <stdlib.h>
#include <math.h>

void s2let_mw_alm2map(complex double* f, const complex double* flm, int L) {
  int spin = 0;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;
  ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);
}

void s2let_mw_map2alm(complex double* flm, const complex double* f, int L) {
  int spin = 0;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;
  ssht_core_mw_forward_sov_conv_sym(flm, f, L, spin, dl_method, verbosity);
}

void s2let_mw_alm2map_real(double* f, const complex double* flm, int L) {
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;
  ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);
}

void s2let_mw_map2alm_real(complex double* flm, const double* f, int L) {
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;
  ssht_core_mw_forward_sov_conv_sym_real(flm, f, L, dl_method, verbosity);
}

void s2let_mw_allocate(complex double **f, int L)
{
  *f = calloc(L * (2*L-1), sizeof **f);
}

void s2let_mw_allocate_real(double **f, int L)
{
  *f = calloc(L * (2*L-1), sizeof **f);
}

double s2let_mw_power(complex double *f, int L){
  complex double *flm;
  s2let_lm_allocate(&flm, L);
  s2let_mw_map2alm(flm, f, L);
  double res = s2let_lm_power(flm, L);
  free(flm);
  return res;
}

double s2let_mw_power_real(double *f, int L){
  complex double *flm;
  s2let_lm_allocate(&flm, L);
  s2let_mw_map2alm_real(flm, f, L);
  double res = s2let_lm_power(flm, L);
  free(flm);
  return res;
}

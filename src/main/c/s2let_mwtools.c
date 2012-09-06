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
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;
  ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);
}

void s2let_mw_map2alm(complex double* flm, const complex double* f, int L) {
  int spin = 0;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;
  ssht_core_mw_forward_sov_conv_sym(flm, f, L, spin, dl_method, verbosity);
}

void s2let_mw_alm2map_real(double* f, const complex double* flm, int L) {
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;
  ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);
}

void s2let_mw_map2alm_real(complex double* flm, const double* f, int L) {
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;
  ssht_core_mw_forward_sov_conv_sym_real(flm, f, L, dl_method, verbosity);
}

void s2let_allocate_lm(complex double **flm, int L)
{
  *flm = (complex double*)calloc(L * L, sizeof(complex double));
}

void s2let_allocate_mw(complex double **f, int L)
{
  *f = (complex double*)calloc(L * (2*L-1), sizeof(complex double));
}

void s2let_allocate_mw_real(double **f, int L)
{
  *f = (double*)calloc(L * (2*L-1), sizeof(double));
}

double s2let_power_lm(complex double *flm, int L){
  int i;
  double totalpower = 0.0;
  for(i = 0; i < L*L; i++)
    totalpower += pow(cabs(flm[i]), 2.0);
  totalpower = totalpower / (L * L);
  return totalpower; 
}

double s2let_power_mw(complex double *f, int L){
  complex double *flm;
  s2let_allocate_lm(&flm, L);
  s2let_mw_map2alm(flm, f, L);
  double res = s2let_power_lm(flm, L);
  free(flm);
  return res; 
}

double s2let_power_mw_real(double *f, int L){
  complex double *flm;
  s2let_allocate_lm(&flm, L);
  s2let_mw_map2alm_real(flm, f, L);
  double res = s2let_power_lm(flm, L);
  free(flm);
  return res; 
}

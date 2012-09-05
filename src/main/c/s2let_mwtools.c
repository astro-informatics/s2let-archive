// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <ssht.h>

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

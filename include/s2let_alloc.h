// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_ALLOC
#define S2LET_ALLOC

#include <complex.h>

/** Pixel space allocation **/

void s2let_mw_allocate(complex double **f, int L);
void s2let_mw_allocate_real(double **f, int L);

void s2let_mwss_allocate(complex double **f, int L);
void s2let_mwss_allocate_real(double **f, int L);

/** Harmonic space allocation **/

void s2let_lm_allocate(complex double **flm, int L);

/** Wigner space allocation **/

void s2let_allocate_f_wav_lmn(
    complex double **f_wav_lmn,
    complex double **f_scal_lm,
    const s2let_parameters_t *parameters
);

/** Wavelet space allocation **/

void s2let_allocate_mw_f_wav(
    complex double **f_wav,
    complex double **f_scal,
    const s2let_parameters_t *parameters
);

void s2let_allocate_mw_f_wav_real(
    double **f_wav,
    double **f_scal,
    const s2let_parameters_t *parameters
);

#endif

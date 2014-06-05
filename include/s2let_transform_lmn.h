// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_TRANSFORM_LMN
#define S2LET_TRANSFORM_LMN

#include <complex.h>

/** Allocation helpers **/

void s2let_allocate_f_wav_lmn(
    complex double **f_wav_lmn,
    complex double **f_scal_lm,
    const s2let_parameters_t *parameters
);

void s2let_allocate_f_wav_lmn_multires(
    complex double **f_wav_lmn,
    complex double **f_scal_lm,
    const s2let_parameters_t *parameters
);

/** Harmonic-space wavelet transform **/

void s2let_wav_analysis_harmonic(
    complex double *f_wav_lmn,
    complex double *f_scal_lm,
    const complex double *flm,
    const complex double *wav_lm,
    const double *scal_l,
    int B,
    int L,
    int J_min,
    int N,
    int spin
);

void s2let_wav_synthesis_harmonic(
    complex double *flm,
    const complex double *f_wav_lmn,
    const complex double *f_scal_lm,
    const complex double *wav_lm,
    const double *scal_l,
    int B,
    int L,
    int J_min,
    int N,
    int spin
);

/** Multi-resolution harmonic-space wavelet transform **/

void s2let_wav_analysis_harmonic_multires(
    complex double *f_wav_lmn,
    complex double *f_scal_lm,
    const complex double *flm,
    const complex double *wav_lm,
    const double *scal_l,
    int B,
    int L,
    int J_min,
    int N,
    int spin
);

void s2let_wav_synthesis_harmonic_multires(
    complex double *flm,
    const complex double *f_wav_lmn,
    const complex double *f_scal_lm,
    const complex double *wav_lm,
    const double *scal_l,
    int B,
    int L,
    int J_min,
    int N,
    int spin
);

#endif

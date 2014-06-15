// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_ALLOC
#define S2LET_ALLOC

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

void s2let_allocate_mw_f_wav(
    complex double **f_wav,
    complex double **f_scal,
    const s2let_parameters_t *parameters
);

void s2let_allocate_mw_f_wav_multires(
    complex double **f_wav,
    complex double **f_scal,
    const s2let_parameters_t *parameters
);

void s2let_allocate_mw_f_wav_real(
    double **f_wav,
    double **f_scal,
    const s2let_parameters_t *parameters
);

void s2let_allocate_mw_f_wav_multires_real(
    double **f_wav,
    double **f_scal,
    const s2let_parameters_t *parameters
);

#endif

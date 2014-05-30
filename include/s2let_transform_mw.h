// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_TRANSFORM_MW
#define S2LET_TRANSFORM_MW

#include <complex.h>

/** Allocation helpers **/

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

/** Pixel-space wavelet transform **/

void s2let_wav_analysis_mw(
    complex double *f_wav,
    complex double *f_scal,
    const complex double *f,
    int B,
    int L,
    int J_min,
    int N,
    int spin,
    s2let_wav_norm_t normalization,
    int original_spin
);

void s2let_wav_synthesis_mw(
    complex double *f,
    const complex double *f_wav,
    const complex double *f_scal,
    int B,
    int L,
    int J_min,
    int N,
    int spin,
    s2let_wav_norm_t normalization,
    int original_spin
);

void s2let_wav_analysis_mw_real(
    double *f_wav,
    double *f_scal,
    const double *f,
    int B,
    int L,
    int J_min,
    int N
);

void s2let_wav_synthesis_mw_real(
    double *f,
    const double *f_wav,
    const double *f_scal,
    int B,
    int L,
    int J_min,
    int N
);

/** Multi-resolution harmonic-space wavelet transform **/

void s2let_wav_analysis_mw_multires(
    complex double *f_wav,
    complex double *f_scal,
    const complex double *f,
    int B,
    int L,
    int J_min,
    int N,
    int spin,
    s2let_wav_norm_t normalization,
    int original_spin
);

void s2let_wav_synthesis_mw_multires(
    complex double *f,
    const complex double *f_wav,
    const complex double *f_scal,
    int B,
    int L,
    int J_min,
    int N,
    int spin,
    s2let_wav_norm_t normalization,
    int original_spin
);

void s2let_wav_analysis_mw_multires_real(
    double *f_wav,
    double *f_scal,
    const double *f,
    int B,
    int L,
    int J_min,
    int N
);

void s2let_wav_synthesis_mw_multires_real(
    double *f,
    const double *f_wav,
    const double *f_scal,
    int B,
    int L,
    int J_min,
    int N
);

// void s2let_axisym_mw_wav_hardthreshold_multires_real(double *g_wav, const double *treshold, int B, int L, int J_min);
// void s2let_axisym_mw_wav_hardthreshold_real(double *g_wav, const double *treshold, int B, int L, int J_min);

#endif

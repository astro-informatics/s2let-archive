// S2LET package
// Copyright (C) 2012 - 2014
// Boris Leistedt & Jason McEwen

#ifndef S2LET_TRANSFORM_LM2WAV
#define S2LET_TRANSFORM_LM2WAV

#include <complex.h>

/** Harmonic-space wavelet transform **/

void s2let_wav_analysis_lm2wav(
    complex double *f_wav,
    complex double *f_scal,
    const complex double *flm,
    const s2let_parameters_t *parameters
);

void s2let_wav_synthesis_lm2wav(
    complex double *flm,
    const complex double *f_wav,
    const complex double *f_scal,
    const s2let_parameters_t *parameters
);

/** Harmonic-space wavelet transform for real signals **/

void s2let_wav_analysis_lm2wav_real(
    double *f_wav,
    double *f_scal,
    const complex double *flm,
    const s2let_parameters_t *parameters
);

void s2let_wav_synthesis_lm2wav_real(
    complex double *flm,
    const double *f_wav,
    const double *f_scal,
    const s2let_parameters_t *parameters
);

/** Multi-resolution harmonic-space wavelet transform **/

void s2let_wav_analysis_lm2wav_multires(
    complex double *f_wav,
    complex double *f_scal,
    const complex double *flm,
    const s2let_parameters_t *parameters
);

void s2let_wav_synthesis_lm2wav_multires(
    complex double *flm,
    const complex double *f_wav,
    const complex double *f_scal,
    const s2let_parameters_t *parameters
);

/** Multi-resolution harmonic-space wavelet transform for real signals **/

void s2let_wav_analysis_lm2wav_multires_real(
    double *f_wav,
    double *f_scal,
    const complex double *flm,
    const s2let_parameters_t *parameters
);

void s2let_wav_synthesis_lm2wav_multires_real(
    complex double *flm,
    const double *f_wav,
    const double *f_scal,
    const s2let_parameters_t *parameters
);

#endif

// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_transform_AXISYM_LM
#define S2LET_transform_AXISYM_LM

#include <complex.h>


void s2let_transform_axisym_lm_allocate_f_wav(
    complex double **f_wav_lm,
    complex double **f_scal_lm,
    const s2let_parameters_t *parameters
);

void s2let_transform_axisym_lm_allocate_f_wav_multires(
    complex double **f_wav_lm,
    complex double **f_scal_lm,
    const s2let_parameters_t *parameters
);

void s2let_transform_axisym_lm_allocate_wav(double **wav_lm, double **scal_lm, const s2let_parameters_t *parameters);
void s2let_transform_axisym_lm_wav(double *wav_lm, double *scal_lm, int B, int L, int J_min);

void s2let_transform_axisym_lm_wav_analysis(
    complex double *f_wav_lm,
    complex double *f_scal_lm,
    const complex double *flm,
    const double *wav_lm,
    const double *scal_lm,
    int B, int L, int J_min
);
void s2let_transform_axisym_lm_wav_synthesis(
    complex double *flm,
    const complex double *f_wav_lm,
    const complex double *f_scal_lm,
    const double *wav_lm,
    const double *scal_lm,
    int B, int L, int J_min
);

void s2let_transform_axisym_lm_wav_analysis_multires(
    complex double *f_wav_lm,
    complex double *f_scal_lm,
    const complex double *flm,
    const double *wav_lm,
    const double *scal_lm,
    int B, int L, int J_min
);
void s2let_transform_axisym_lm_wav_synthesis_multires(
    complex double *flm,
    const complex double *f_wav_lm,
    const complex double *f_scal_lm,
    const double *wav_lm,
    const double *scal_lm,
    int B, int L, int J_min
);

#endif

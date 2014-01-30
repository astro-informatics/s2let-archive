// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_TRANSFORM_LMN
#define S2LET_TRANSFORM_LMN

#include <complex.h>

void s2let_allocate_f_wav_lmn(
    complex double **f_wav_lmn,
    complex double **f_scal_lm,
    int B,
    int L,
    int J_min,
    int N
);

/*void s2let_allocate_f_wav_lmn_multires(
    complex double **f_wav_lmn,
    complex double **f_scal_lmn,
    int B,
    int L,
    int J_min,
    int N
);*/


void s2let_wav_analysis_harmonic(
    complex double *f_wav_lmn,
    complex double *f_scal_lm,
    const complex double *flm,
    const complex double *wav_lm,
    const double *scal_l,
    int B,
    int L,
    int J_min,
    int N
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
    int N
);

//void s2let_axisym_lm_wav_analysis_multires(complex double *f_wav_lm, complex double *f_scal_lm, const complex double *flm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min);
//void s2let_axisym_lm_wav_synthesis_multires(complex double *flm, const complex double *f_wav_lm, const complex double *f_scal_lm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min);


#endif

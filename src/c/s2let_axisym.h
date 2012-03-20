// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_AXISYM
#define S2LET_AXISYM

void s2let_axisym_allocate_wav_lm(double **wav_lm, double **scal_lm, int B, int L);
void s2let_axisym_wav_lm(double *wav_lm, double *scal_lm, int B, int L, int J_min);

void s2let_axisym_allocate_f_wav_lm(complex double **f_wav_lm, complex double **f_scal_lm, int B, int L);
void s2let_axisym_allocate_f_wav_lm_real(double **f_wav_lm, double **f_scal_lm, int B, int L);

void s2let_axisym_wav_analysis_lm(complex double *f_wav_lm, complex double *f_scal_lm, const complex double *flm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min);
void s2let_axisym_wav_synthesis_lm(complex double *flm, const complex double *f_wav_lm, const complex double *f_scal_lm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min);

void s2let_axisym_wav_analysis(complex double *f_wav, complex double *f_scal, const complex double *f, int B, int L, int J_min);
void s2let_axisym_wav_analysis_real(double *f_wav, double *f_scal, const double *f, int B, int L, int J_min);
void s2let_axisym_wav_synthesis(complex double *f, const complex double *f_wav, const complex double *f_scal, int B, int L, int J_min);
void s2let_axisym_wav_synthesis_real(double *f, const double *f_wav, const double *f_scal, int B, int L, int J_min);

#endif
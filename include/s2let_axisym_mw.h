// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_AXISYM_MW
#define S2LET_AXISYM_MW

#include <complex.h>

void s2let_axisym_mw_allocate_f_wav_multires(complex double **f_wav, complex double **f_scal, int B, int L, int J_min);
void s2let_axisym_mw_allocate_f_wav_multires_real(double **f_wav, double **f_scal, int B, int L, int J_min);
void s2let_axisym_mw_allocate_f_wav(complex double **f_wav, complex double **f_scal, int B, int L, int J_min);
void s2let_axisym_mw_allocate_f_wav_real(double **f_wav, double **f_scal, int B, int L, int J_min);

void s2let_axisym_mw_wav_analysis(complex double *f_wav, complex double *f_scal, const complex double *f, int B, int L, int J_min);
void s2let_axisym_mw_wav_synthesis(complex double *f, const complex double *f_wav, const complex double *f_scal, int B, int L, int J_min);

void s2let_axisym_mw_wav_analysis_real(double *f_wav, double *f_scal, const double *f, int B, int L, int J_min);
void s2let_axisym_mw_wav_synthesis_real(double *f, const double *f_wav, const double *f_scal, int B, int L, int J_min);

void s2let_axisym_mw_wav_analysis_multires(complex double *f_wav, complex double *f_scal, const complex double *f, int B, int L, int J_min);
void s2let_axisym_mw_wav_synthesis_multires(complex double *f, const complex double *f_wav, const complex double *f_scal, int B, int L, int J_min);

void s2let_axisym_mw_wav_analysis_multires_real(double *f_wav, double *f_scal, const double *f, int B, int L, int J_min);
void s2let_axisym_mw_wav_synthesis_multires_real(double *f, const double *f_wav, const double *f_scal, int B, int L, int J_min);

void s2let_axisym_mw_wav_hardthreshold_multires_real(double *g_wav, const double *threshold, int B, int L, int J_min);
void s2let_axisym_mw_wav_hardthreshold_real(double *g_wav, const double *threshold, int B, int L, int J_min);

#endif

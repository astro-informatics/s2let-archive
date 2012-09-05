// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_AXISYM_HPX
#define S2LET_AXISYM_HPX

void s2let_axisym_hpx_allocate_f_wav_real(double **f_wav, double **f_scal, int nside, int B, int L, int J_min);

void s2let_axisym_hpx_wav_analysis_real(double *f_wav, double *f_scal, const double *f, int nside, int B, int L, int J_min);
void s2let_axisym_hpx_wav_synthesis_real(double *f, const double *f_wav, const double *f_scal, int nside, int B, int L, int J_min);

#endif
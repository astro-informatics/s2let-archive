// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_CORE
#define S2LET_CORE

#define PI    3.141592653589793238462643383279502884197

int s2let_el_min(int B, int J_min);

int s2let_j_max(int L, int B);

void s2let_Lj(int *Lj, int B, int J_min, int L);

void allocate_tilling(double **kappa, double **kappa0, int B, int L);

void s2let_tilling(double *kappa, double *kappa0, int B, int L);

void s2let_allocate_wav_lm(double **wav_lm, double **scal_lm, int B, int L);

void s2let_wav_lm(double *wav_lm, double *scal_lm, int B, int L, int J_min);

#endif
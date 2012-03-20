// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_TILLING
#define S2LET_TILLING

#define PI    3.141592653589793238462643383279502884197

int s2let_el_min(int B, int J_min);

int s2let_j_max(int L, int B);

void s2let_Lj(int *Lj, int B, int L, int J_min);

void s2let_allocate_tilling(double **kappa, double **kappa0, int B, int L);

void s2let_tilling(double *kappa, double *kappa0, int B, int L, int J_min);

void s2let_tilling_phi2(double *phi2, int B, int L, int J_min);

double s2let_check_identity(double *kappa, double *kappa0, int B, int L, int J_min);

int jlm2ind(int j, int l, int m, int L);

int lm2ind(int l, int m);

#endif
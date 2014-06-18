// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_LM
#define S2LET_LM

double s2let_lm_power(complex double *flm, int L);

void s2let_lm_random_flm(complex double *flm, int L, int spin, int seed);
void s2let_lm_random_flm_real(complex double *flm, int L, int seed);


#endif

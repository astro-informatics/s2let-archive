// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_MATH
#define S2LET_MATH

#define PI    3.141592653589793238462643383279502884197

double s2let_math_kappa0_quadtrap_s2dw(double a, double b, int n, int B);

double s2let_math_kappa0_quadtrap_needlet(double a, double b, int n);

/*!
 * Random number from seed (Numerical Recipes)
 */
double ran2_dp(int idum);

/*!
 * Max absolute error between two complex arrays
 */
double maxerr_cplx(complex double *a, complex double *b, int size);

/*!
 * Max absolute error between two real arrays
 */
double maxerr(double *a, double *b, int size);


/*!
 * Computes power of a signal from its harmonics
 */
double s2let_power_lm(complex double *flm, int L);

#endif

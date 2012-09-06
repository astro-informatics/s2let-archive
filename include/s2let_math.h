// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_MATH
#define S2LET_MATH

#define PI    3.141592653589793238462643383279502884197

double s2let_kappa0_integrand(double k, int B);

double s2let_kappa0_quadtrap(double a, double b, int n, int B);

/*!
 * Random number from seed (Numerical Recipes)
 */
double ran2_dp(int idum);

/*!
 * Generate random harmonic coefficients for a complex map.
 *
 * \param[out]  flm Harmonic coefficients.
 * \param[in]  L Band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_axisym_random_flm(complex double *flm, int L, int seed);

/*!
 * Generate random harmonic coefficients corresponding to a real map.
 *
 * \param[out]  flm Harmonic coefficients.
 * \param[in]  L Band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_axisym_random_flm_real(complex double *flm, int L, int seed);

/*!
 * Max absolute error between two complex arrays
 */
double maxerr_cplx(complex double *a, complex double *b, int size);

/*!
 * Max absolute error between two real arrays
 */
double maxerr(double *a, double *b, int size);

#endif

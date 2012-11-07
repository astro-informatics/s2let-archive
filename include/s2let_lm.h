// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_LM
#define S2LET_LM


double s2let_lm_power(complex double *flm, int L);


/*!
 * Generate random harmonic coefficients for a complex map.
 *
 * \param[out]  flm Harmonic coefficients.
 * \param[in]  L Band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_lm_random_flm(complex double *flm, int L, int seed);

/*!
 * Generate random harmonic coefficients corresponding to a real map.
 *
 * \param[out]  flm Harmonic coefficients.
 * \param[in]  L Band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_lm_random_flm_real(complex double *flm, int L, int seed);


/*!
 * Allocate spherical har for a given bandlimit L
 *
 * \param[inout]  flm Spherical harmonic coefficients.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_lm_allocate(complex double **flm, int L);


#endif

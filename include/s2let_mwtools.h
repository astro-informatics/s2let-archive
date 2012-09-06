// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_MWTOOLS
#define S2LET_MWTOOLS

/*!
 * Interface to SSHT (required by the Java interface to S2LET)
 * Output : MW complex signal
 *
 * \param[out]  f Output MW map.
 * \param[in]  flm Spherical harmonic coefficients.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_mw_alm2map(complex double* f, const complex double* flm, int L);

/*!
 * Interface to SSHT (required by the Java interface to S2LET)
 * Output : spherical harmonics
 *
 * \param[out]  flm Spherical harmonic coefficients.
 * \param[in]  f Output MW map.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_mw_map2alm(complex double* flm, const complex double* f, int L);

/*!
 * Interface to SSHT (required by the Java interface to S2LET)
 * Output : MW real signal
 *
 * \param[out]  f Output MW map.
 * \param[in]  flm Spherical harmonic coefficients.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_mw_alm2map_real(double* f, const complex double* flm, int L);

/*!
 * Interface to SSHT (required by the Java interface to S2LET)
 * Output : spherical harmonics (corresponding to real signal)
 *
 * \param[out]  flm Spherical harmonic coefficients.
 * \param[in]  f Output MW map.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_mw_map2alm_real(complex double* flm, const double* f, int L);

/*!
 * Allocate spherical har for a given bandlimit L
 *
 * \param[inout]  flm Spherical harmonic coefficients.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_allocate_lm(complex double **flm, int L);

/*!
 * Allocate MW map for a given bandlimit L
 *
 * \param[inout]  f MW map to be allocated.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_allocate_mw(complex double **f, int L);

/*!
 * Allocate real MW map for a given bandlimit L
 *
 * \param[inout]  f Real MW map to be allocated.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_allocate_mw_real(double **f, int L);

/*!
 * Computes power of a signal from its harmonics
 */
double s2let_power_lm(complex double *flm, int L);

/*!
 * Computes power of a MW signal
 */
double s2let_power_mw(complex double *flm, int L);

/*!
 * Computes power of a MW signal
 */
double s2let_power_mw_real(double *flm, int L);

#endif

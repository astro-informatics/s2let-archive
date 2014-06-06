// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <complex.h>
#include <ssht.h>
#include <so3.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


double s2let_lm_power(complex double *flm, int L){
    int i;
    double totalpower = 0.0;
    for(i = 0; i < L*L; ++i)
        totalpower += pow(cabs(flm[i]), 2.0);
    totalpower = totalpower / (L * L);
    return totalpower;
}

/*!
 * Allocate spherical harmonic coefficients for a given
 * bandlimit L.
 *
 * \param[out]  flm Pointer to allocated space for spherical
 *                  harmonic coefficients.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_lm_allocate(complex double **flm, int L)
{
    *flm = calloc(L * L, sizeof **flm);
}

/*!
 * Allocate Wigner coefficients for given bandlimits L and N.
 *
 * \param[out]  flmn Pointer to allocated space for Wigner
 *                   coefficients.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Azimuthal harmonic band-limit.
 * \retval none
 */
void s2let_lmn_allocate(complex double **flmn, const so3_parameters_t *parameters)
{
    // Create a copy of the parameters where reality is set to false, in order to
    // always allocate the full lmn block.
    // TODO: Make sure that lmn and mw transforms can actually deal will the
    // reduced storage for real signals, and then get rid of this copy here.
    so3_parameters_t complex_parameters = *parameters;
    complex_parameters.reality = 0;

    *flmn = calloc(so3_sampling_flmn_size(&complex_parameters), sizeof **flmn);
}

/*!
 * Generate random harmonic coefficients for a complex map.
 *
 * \param[out]  flm Harmonic coefficients.
 * \param[in]  L Band-limit.
 * \parma[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_lm_random_flm(complex double *flm, int L, int spin, int seed)
{
    int i, i_min;
    srand( time(NULL) );
    // el < |s| are zero, so start with el = |s|, m = -el.
    i_min = spin*spin;
    for (i=i_min; i<L*L; ++i)
        flm[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
}

/*!
 * Generate random harmonic coefficients corresponding to a real map.
 *
 * \param[out]  flm Harmonic coefficients.
 * \param[in]  L Band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_lm_random_flm_real(complex double *flm, int L, int seed) {
    int el, m, msign, i, i_op;
    for (el = 0; el < L; ++el) {
        m = 0;
        i = el*el + el + m;
        flm[i] = (2.0*ran2_dp(seed) - 1.0);
        for (m = 1; m <= el; ++m) {
            i = el*el + el + m;
            flm[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
            i_op = el*el + el - m;
            msign = m & 1;
            msign = 1 - msign - msign; // (-1)^m
            flm[i_op] = msign * conj(flm[i]);
        }
    }
}

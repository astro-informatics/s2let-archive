// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <so3.h>
#include <ssht.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

static inline void fill_so3_parameters(so3_parameters_t *parameters, int L, int N)
{
    parameters->L = L;
    parameters->N = N;
    parameters->n_order = SO3_N_ORDER_NEGATIVE_FIRST;
    parameters->storage = SO3_STORAGE_PADDED;
}

/*!
 * Allocates arrays for directional wavelet transform in harmonic space.
 *
 * \param[out]  f_wav_lmn Wigner coefficients of the wavelet contributions.
 *                        Each wavelet has size (2*N-1)*L*L and there are
 *                        (J-J_min+1) scales.
 * \param[out]  f_scal_lm Spherical harmonic coefficients of the scaling
 *                        contribution (L*L).
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \retval none
 */
void s2let_allocate_f_wav_lmn(
    complex double **f_wav_lmn,
    complex double **f_scal_lm,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;

    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, L, N);

    int lmn_block_size = so3_sampling_flmn_size(&so3_parameters);

    int J = s2let_j_max(parameters);
    *f_wav_lmn = calloc((J-J_min+1) * lmn_block_size, sizeof **f_wav_lmn);
    *f_scal_lm = calloc(L*L, sizeof **f_scal_lm);
}

/*!
 * Allocates multi-resolution arrays for directional wavelet transform in
 * harmonic space.
 *
 * \param[out]  f_wav_lmn Wigner coefficients of the wavelet contributions.
 *                        The size of each wavelet depends on its band-limit.
 * \param[out]  f_scal_lm Spherical harmonic coefficients of the scaling
 *                        contribution (L*L).
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \retval none
 */
void s2let_allocate_f_wav_lmn_multires(
    complex double **f_wav_lmn,
    complex double **f_scal_lm,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;

    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, L, N);

    int J = s2let_j_max(parameters);
    int j, bandlimit, total = 0;
    for (j = J_min; j <= J; ++j)
    {
        bandlimit = MIN(s2let_bandlimit(j, parameters), L);
        so3_parameters.L = bandlimit;
        total += so3_sampling_flmn_size(&so3_parameters);
    }
    *f_wav_lmn = calloc(total, sizeof **f_wav_lmn);
    bandlimit = MIN(s2let_bandlimit(J_min-1, parameters), L);
    *f_scal_lm = calloc(bandlimit * bandlimit, sizeof **f_scal_lm);
}

/*!
 * Spherical wavelets: full resolution analysis in harmonic space.
 * Perform directional wavelet transform from precomputed kernels
 * to give Wigner coefficients.
 *
 * \param[out]  f_wav_lmn Wavelet transform (Wigner coefficients of wavelet contribution).
 * \param[out]  f_scal_lm Wavelet transform (spherical harmonic coefficients of scaling contribution).
 * \param[in]  flm Spherical harmonic decomposition of input function.
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_l Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \retval none
 */
void s2let_wav_analysis_harmonic(
    complex double *f_wav_lmn,
    complex double *f_scal_lm,
    const complex double *flm,
    const complex double *wav_lm,
    const double *scal_l,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;
    int spin = parameters->spin;

    int j, el, m ,n;
    int J = s2let_j_max(parameters);

    int lm_ind, lmn_ind;
    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, L, N);

    complex double psi;
    double phi;

    int offset = 0;

    for (j = J_min; j <= J; ++j)
    {
        for (n = -N+1; n < N; ++n)
        {
            for (el = MAX(ABS(spin), ABS(n)); el < L; ++el)
            {
                ssht_sampling_elm2ind(&lm_ind, el, n);
                psi = 8*PI*PI/(2*el+1) * conj(wav_lm[j*L*L + lm_ind]);
                for (m = -el; m <= el; ++m)
                {
                    ssht_sampling_elm2ind(&lm_ind, el, m);
                    so3_sampling_elmn2ind(&lmn_ind, el, m, n, &so3_parameters);
                    f_wav_lmn[offset + lmn_ind] = flm[lm_ind] * psi;
                }
            }
        }
        offset += so3_sampling_flmn_size(&so3_parameters);
    }

    for (el = ABS(spin); el < L; ++el)
    {
        phi = sqrt(4.0*PI/(2*el+1)) * scal_l[el];
        for (m = -el; m <= el; ++m)
        {
            ssht_sampling_elm2ind(&lm_ind, el, m);
            f_scal_lm[lm_ind] = flm[lm_ind] * phi;
        }
    }
}


/*!
 * Spherical wavelets: full resolution synthesis in harmonic space.
 * Perform directional wavelet transform from precomputed kernels
 * to give Wigner coefficients.
 *
 * \param[out]  flm Spherical harmonic decomposition of input function.
 * \param[in]  f_wav_lmn Wavelet transform (Wigner coefficients of wavelet contribution).
 * \param[in]  f_scal_lm Wavelet transform (spherical harmonic coefficients of scaling contribution).
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_l Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \retval none
 */
void s2let_wav_synthesis_harmonic(
    complex double *flm,
    const complex double *f_wav_lmn,
    const complex double *f_scal_lm,
    const complex double *wav_lm,
    const double *scal_l,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;
    int spin = parameters->spin;

    int j, el, m ,n;
    int J = s2let_j_max(parameters);

    int lm_ind, lmn_ind;
    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, L, N);

    complex double psi;
    double phi;

    int offset = 0;


    for (j = J_min; j <= J; ++j)
    {
        for (n = -N+1; n < N; ++n)
        {
            for (el = MAX(ABS(spin), ABS(n)); el < L; ++el)
            {
                ssht_sampling_elm2ind(&lm_ind, el, n);
                psi = wav_lm[j*L*L + lm_ind];
                for (m = -el; m <= el; ++m)
                {
                    ssht_sampling_elm2ind(&lm_ind, el, m);
                    so3_sampling_elmn2ind(&lmn_ind, el, m, n, &so3_parameters);
                    flm[lm_ind] += f_wav_lmn[offset + lmn_ind] * psi;
                }
            }
        }
        offset += so3_sampling_flmn_size(&so3_parameters);
    }

    for (el = ABS(spin); el < L; ++el)
    {
        phi = sqrt(4.0*PI/(2*el+1)) * scal_l[el];
        for (m = -el; m <= el; ++m)
        {
            ssht_sampling_elm2ind(&lm_ind, el, m);
            flm[lm_ind] += f_scal_lm[lm_ind] * phi;
        }
    }
}

/*!
 * Spherical wavelets: multi-resolution analysis in harmonic space.
 * Perform directional wavelet transform from precomputed kernels
 * to give Wigner coefficients.
 *
 * \param[out]  f_wav_lmn Wavelet transform (Wigner coefficients of wavelet contribution).
 * \param[out]  f_scal_lm Wavelet transform (spherical harmonic coefficients of scaling contribution).
 * \param[in]  flm Spherical harmonic decomposition of input function.
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_l Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \retval none
 */
void s2let_wav_analysis_harmonic_multires(
    complex double *f_wav_lmn,
    complex double *f_scal_lm,
    const complex double *flm,
    const complex double *wav_lm,
    const double *scal_l,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;
    int spin = parameters->spin;

    int j, el, m ,n;
    int J = s2let_j_max(parameters);
    int bandlimit;

    int lm_ind, lmn_ind;
    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, L, N);

    complex double psi;
    double phi;

    int offset = 0;

    for (j = J_min; j <= J; j++)
    {
        bandlimit = MIN(s2let_bandlimit(j, parameters), L);
        so3_parameters.L = bandlimit;

        for (n = -N+1; n < N; ++n)
        {
            for (el = MAX(ABS(spin), ABS(n)); el < bandlimit; ++el)
            {
                ssht_sampling_elm2ind(&lm_ind, el, n);
                psi = 8*PI*PI/(2*el+1) * conj(wav_lm[j*L*L + lm_ind]);
                for (m = -el; m <= el; ++m)
                {
                    ssht_sampling_elm2ind(&lm_ind, el, m);
                    so3_sampling_elmn2ind(&lmn_ind, el, m, n, &so3_parameters);
                    f_wav_lmn[offset + lmn_ind] = flm[lm_ind] * psi;
                }
            }
        }
        offset += so3_sampling_flmn_size(&so3_parameters);
    }

    bandlimit = MIN(s2let_bandlimit(J_min-1, parameters), L);
    for (el = ABS(spin); el < bandlimit; ++el)
    {
        phi = sqrt(4.0*PI/(2*el+1)) * scal_l[el];
        for (m = -el; m <= el; ++m)
        {
            ssht_sampling_elm2ind(&lm_ind, el, m);
            f_scal_lm[lm_ind] = flm[lm_ind] * phi;
        }
    }
}


/*!
 * Spherical wavelets: multi-resolution synthesis in harmonic space.
 * Perform directional wavelet transform from precomputed kernels
 * to give Wigner coefficients.
 *
 * \param[out]  flm Spherical harmonic decomposition of input function.
 * \param[in]  f_wav_lmn Wavelet transform (Wigner coefficients of wavelet contribution).
 * \param[in]  f_scal_lm Wavelet transform (spherical harmonic coefficients of scaling contribution).
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_l Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \retval none
 */
void s2let_wav_synthesis_harmonic_multires(
    complex double *flm,
    const complex double *f_wav_lmn,
    const complex double *f_scal_lm,
    const complex double *wav_lm,
    const double *scal_l,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;
    int spin = parameters->spin;

    int j, el, m ,n;
    int J = s2let_j_max(parameters);
    int bandlimit;

    int lm_ind, lmn_ind;
    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, L, N);

    complex double psi;
    double phi;

    int offset = 0;

    for (j = J_min; j <= J; ++j)
    {
        bandlimit = MIN(s2let_bandlimit(j, parameters), L);
        so3_parameters.L = bandlimit;

        for (n = -N+1; n < N; ++n)
        {
            for (el = MAX(ABS(spin), ABS(n)); el < bandlimit; ++el)
            {
                ssht_sampling_elm2ind(&lm_ind, el, n);
                psi = wav_lm[j*L*L + lm_ind];
                for (m = -el; m <= el; ++m)
                {
                    ssht_sampling_elm2ind(&lm_ind, el, m);
                    so3_sampling_elmn2ind(&lmn_ind, el, m, n, &so3_parameters);
                    flm[lm_ind] += f_wav_lmn[offset + lmn_ind] * psi;
                }
            }
        }
        offset += so3_sampling_flmn_size(&so3_parameters);
    }

    bandlimit = MIN(s2let_bandlimit(J_min-1, parameters), L);
    for (el = ABS(spin); el < bandlimit; ++el)
    {
        phi = sqrt(4.0*PI/(2*el+1)) * scal_l[el];
        for (m = -el; m <= el; ++m)
        {
            ssht_sampling_elm2ind(&lm_ind, el, m);
            flm[lm_ind] += f_scal_lm[lm_ind] * phi;
        }
    }
}

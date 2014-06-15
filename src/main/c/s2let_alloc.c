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

static inline void fill_so3_parameters(so3_parameters_t *so3_parameters, const s2let_parameters_t *parameters)
{
    so3_parameters->verbosity = parameters->verbosity;
    so3_parameters->L = parameters->L;
    so3_parameters->N = parameters->N;
    so3_parameters->sampling_scheme = parameters->sampling_scheme;
    so3_parameters->n_order = S2LET_SO3_N_ORDER;
    so3_parameters->storage = S2LET_SO3_STORAGE;
    so3_parameters->dl_method = parameters->dl_method;
    so3_parameters->reality = parameters->reality;

    if (parameters->N % 2)
        so3_parameters->n_mode = SO3_N_MODE_EVEN;
    else
        so3_parameters->n_mode = SO3_N_MODE_ODD;
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

    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, parameters);

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
    fill_so3_parameters(&so3_parameters, parameters);

    int J = s2let_j_max(parameters);
    int j, bandlimit, total = 0;
    for (j = J_min; j <= J; ++j)
    {
        bandlimit = MIN(s2let_bandlimit(j, parameters), L);
        so3_parameters.L = bandlimit;
        so3_parameters.N = MIN(N,bandlimit);
        total += so3_sampling_flmn_size(&so3_parameters);
    }

    *f_wav_lmn = calloc(total, sizeof **f_wav_lmn);
    bandlimit = MIN(s2let_bandlimit(J_min-1, parameters), L);
    *f_scal_lm = calloc(bandlimit * bandlimit, sizeof **f_scal_lm);
}


/*!
 * Allocates arrays for wavelet transform in pixel space (MW sampling).
 *
 * \param[out]  f_wav Pointer to allocated space for array of wavelet
 *                    maps, using MW sampling.
 * \param[out]  f_scal Pointer to allocated space for scaling function
 *                     map, using MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \retval none
 */
void s2let_allocate_mw_f_wav(
    complex double **f_wav,
    complex double **f_scal,
    const s2let_parameters_t *parameters
) {
    int J_min = parameters->J_min;

    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, parameters);

    int f_block_size = so3_sampling_f_size(&so3_parameters);

    int J = s2let_j_max(parameters);
    // We actually only need N samples of the orientational angle.
    *f_wav = calloc((J-J_min+1) * f_block_size, sizeof **f_wav);
    *f_scal = calloc(f_block_size, sizeof **f_scal);
}

/*!
 * Allocates multi-resolution arrays for wavelet transform in pixel space
 * (MW sampling).
 *
 * \param[out]  f_wav Pointer to allocated space for array of wavelet
 *                    maps, using MW sampling.
 * \param[out]  f_scal Pointer to allocated space for scaling function
 *                     map, using MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \retval none
 */
void s2let_allocate_mw_f_wav_multires(
    complex double **f_wav,
    complex double **f_scal,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;

    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, parameters);

    int J = s2let_j_max(parameters);
    int j, bandlimit, total = 0;

    for (j = J_min; j <= J; ++j)
    {
        bandlimit = MIN(s2let_bandlimit(j, parameters), L);
        so3_parameters.L = bandlimit;
        // We actually only need N samples of the orientational angle.
        total += so3_sampling_f_size(&so3_parameters);
    }

    *f_wav = calloc(total, sizeof **f_wav);
    bandlimit = MIN(s2let_bandlimit(J_min-1, parameters), L);
    so3_parameters.L = bandlimit;
    total = so3_sampling_f_size(&so3_parameters);
    *f_scal = calloc(total, sizeof **f_scal);
}

/*!
 * Allocates arrays for wavelet transform of real signal in pixel
 * space (MW sampling).
 *
 * \param[out]  f_wav Pointer to allocated space for array of wavelet
 *                    maps, using MW sampling.
 * \param[out]  f_scal Pointer to allocated space for scaling function
 *                     map, using MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \retval none
 */
void s2let_allocate_mw_f_wav_real(
    double **f_wav,
    double **f_scal,
    const s2let_parameters_t *parameters
) {
    int J_min = parameters->J_min;

    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, parameters);

    int f_block_size = so3_sampling_f_size(&so3_parameters);

    int J = s2let_j_max(parameters);
    // We actually only need N samples of the orientational angle.
    *f_wav = calloc((J-J_min+1) * f_block_size, sizeof **f_wav);
    *f_scal = calloc(f_block_size, sizeof **f_scal);
}

/*!
 * Allocates multi-resolution arrays for wavelet transform of real
 * signal in pixel space (MW sampling).
 *
 * \param[out]  f_wav Pointer to allocated space for array of wavelet
 *                    maps, using MW sampling.
 * \param[out]  f_scal Pointer to allocated space for scaling function
 *                     map, using MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \retval none
 */
void s2let_allocate_mw_f_wav_multires_real(
    double **f_wav,
    double **f_scal,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;

    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, parameters);

    int J = s2let_j_max(parameters);
    int j, bandlimit, total = 0;

    for (j = J_min; j <= J; ++j)
    {
        bandlimit = MIN(s2let_bandlimit(j, parameters), L);
        so3_parameters.L = bandlimit;
        // We actually only need N samples of the orientational angle.
        total += so3_sampling_f_size(&so3_parameters);
    }

    *f_wav = calloc(total, sizeof **f_wav);
    bandlimit = MIN(s2let_bandlimit(J_min-1, parameters), L);
    so3_parameters.L = bandlimit;
    total = so3_sampling_f_size(&so3_parameters);
    *f_scal = calloc(total, sizeof **f_scal);
}

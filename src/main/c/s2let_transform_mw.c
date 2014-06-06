// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <complex.h>
#include <ssht.h>
#include <so3.h>
#include <stdlib.h>

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
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;

    int J = s2let_j_max(parameters);
    // We actually only need N samples of the orientational angle.
    *f_wav = calloc((J-J_min+1) * (2*N-1) * L * (2*L-1), sizeof **f_wav);
    *f_scal = calloc(L * (2*L-1), sizeof **f_scal);
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
    int N = parameters->N;

    int J = s2let_j_max(parameters);
    int j, bandlimit, total = 0;

    for (j = J_min; j <= J; ++j)
    {
        bandlimit = MIN(s2let_bandlimit(j, parameters), L);
        // We actually only need N samples of the orientational angle.
        total += (2*N-1) * bandlimit * (2 * bandlimit - 1);
    }

    *f_wav = calloc(total, sizeof **f_wav);
    bandlimit = MIN(s2let_bandlimit(J_min-1, parameters), L);
    *f_scal = calloc(bandlimit * (2*bandlimit-1), sizeof **f_scal);
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
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;

    int J = s2let_j_max(parameters);
    // We actually only need N samples of the orientational angle.
    *f_wav = calloc((J-J_min+1) * (2*N-1) * L * (2*L-1), sizeof **f_wav);
    *f_scal = calloc(L * (2*L-1), sizeof **f_scal);
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
    int N = parameters->N;

    int J = s2let_j_max(parameters);
    int j, bandlimit, total = 0;

    for (j = J_min; j <= J; ++j)
    {
        bandlimit = MIN(s2let_bandlimit(j, parameters), L);
        // We actually only need N samples of the orientational angle.
        total += (2*N-1) * bandlimit * (2 * bandlimit - 1);
    }

    *f_wav = calloc(total, sizeof **f_wav);
    bandlimit = MIN(s2let_bandlimit(J_min-1, parameters), L);
    *f_scal = calloc(bandlimit * (2*bandlimit-1), sizeof **f_scal);
}

/*!
 * Spherical wavelets : full resolution analysis in real space, MW sampling.
 * Perform directional wavelet transform in real space (from scratch,
 * gives wavelet maps).
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f_wav Array of wavelet maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  f Input function (MW sampling)
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  normalization Indicates how to normalise the wavelets
 *                           and scaling function.
 * \param[in]  original_spin If normalization has value
 *                           S2LET_WAV_NORM_SPIN_LOWERED, this parameter
 *                           indicates which spin number the wavelets
 *                           were lowered from. Otherwise, it is ignored.
 * \retval none
 */
void s2let_wav_analysis_mw(
    complex double *f_wav,
    complex double *f_scal,
    const complex double *f,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;
    int spin = parameters->spin;

    int verbosity = 0;
    ssht_dl_method_t dl_method = SSHT_DL_RISBO;
    so3_parameters_t so3_parameters = {};
    so3_parameters.verbosity = verbosity;
    so3_parameters.L = L;
    so3_parameters.N = N;
    so3_parameters.n_order = SO3_N_ORDER_NEGATIVE_FIRST;
    so3_parameters.storage = SO3_STORAGE_PADDED;
    so3_parameters.n_mode = SO3_N_MODE_ALL;
    so3_parameters.dl_method = dl_method;

    int j, offset, offset_lmn;
    int J = s2let_j_max(parameters);
    //int l_min = s2let_axisym_el_min(B, J_min);

    complex double *wav_lm;
    double *scal_l;
    s2let_tiling_wavelet_allocate(&wav_lm, &scal_l, parameters);
    s2let_tiling_wavelet(wav_lm, scal_l, parameters);

    complex double *flm, *f_wav_lmn, *f_scal_lm;

    s2let_lm_allocate(&flm, L);
    ssht_core_mw_forward_sov_conv_sym(flm, f, L, spin, dl_method, verbosity);

    s2let_allocate_f_wav_lmn(&f_wav_lmn, &f_scal_lm, parameters);
    s2let_wav_analysis_harmonic(f_wav_lmn, f_scal_lm, flm, wav_lm, scal_l, parameters);

    // Note, this is a spin-0 transform!
    ssht_core_mw_inverse_sov_sym(f_scal, f_scal_lm, L, 0, dl_method, verbosity);
    offset = 0;
    offset_lmn = 0;
    for (j = J_min; j <= J; ++j)
    {
        so3_core_inverse_via_ssht(
            f_wav + offset,
            f_wav_lmn + offset_lmn,
            &so3_parameters
        );
        offset_lmn += (2*N-1) * L*L;
        offset += (2*N-1) * L * (2*L-1);
    }

    free(flm);
    free(wav_lm);
    free(scal_l);
    free(f_scal_lm);
    free(f_wav_lmn);
}

/*!
 * Spherical wavelets : full resolution synthesis in real space, MW sampling.
 * Perform directional wavelet transform in real space (from scratch,
 * gives wavelet maps).
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Array of wavelets maps, MW sampling.
 * \param[in]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  normalization Indicates how to normalise the wavelets
 *                           and scaling function.
 * \param[in]  original_spin If normalization has value
 *                           S2LET_WAV_NORM_SPIN_LOWERED, this parameter
 *                           indicates which spin number the wavelets
 *                           were lowered from. Otherwise, it is ignored.
 * \retval none
 */
void s2let_wav_synthesis_mw(
    complex double *f,
    const complex double *f_wav,
    const complex double *f_scal,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;
    int spin = parameters->spin;

    int verbosity = 0;
    ssht_dl_method_t dl_method = SSHT_DL_RISBO;
    so3_parameters_t so3_parameters = {};
    so3_parameters.verbosity = verbosity;
    so3_parameters.L = L;
    so3_parameters.N = N;
    so3_parameters.n_order = SO3_N_ORDER_NEGATIVE_FIRST;
    so3_parameters.storage = SO3_STORAGE_PADDED;
    so3_parameters.n_mode = SO3_N_MODE_ALL;
    so3_parameters.dl_method = dl_method;

    int j, offset, offset_lmn;
    int J = s2let_j_max(parameters);
    //int l_min = s2let_axisym_el_min(B, J_min);

    complex double *wav_lm;
    double *scal_l;
    s2let_tiling_wavelet_allocate(&wav_lm, &scal_l, parameters);
    s2let_tiling_wavelet(wav_lm, scal_l, parameters);

    complex double *flm, *f_wav_lmn, *f_scal_lm;
    s2let_allocate_f_wav_lmn(&f_wav_lmn, &f_scal_lm, parameters);

    // Note, this is a spin-0 transform!
    ssht_core_mw_forward_sov_conv_sym(f_scal_lm, f_scal, L, 0, dl_method, verbosity);
    offset = 0;
    offset_lmn = 0;
    for (j = J_min; j <= J; ++j)
    {
        so3_core_forward_via_ssht(
            f_wav_lmn + offset_lmn,
            f_wav + offset,
            &so3_parameters
        );
        offset_lmn += (2*N-1) * L*L;
        offset += (2*N-1) * L * (2*L-1);
    }

    s2let_lm_allocate(&flm, L);
    s2let_wav_synthesis_harmonic(flm, f_wav_lmn, f_scal_lm, wav_lm, scal_l, parameters);

    ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

    free(flm);
    free(wav_lm);
    free(scal_l);
    free(f_scal_lm);
    free(f_wav_lmn);
}

/*!
 * Spherical wavelets : full resolution analysis in real space, MW sampling.
 * Perform directional wavelet transform in real space (from scratch,
 * gives wavelet maps).
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f_wav Array of wavelet maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  f Input function (MW sampling)
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  normalization Indicates how to normalise the wavelets
 *                           and scaling function.
 * \param[in]  original_spin If normalization has value
 *                           S2LET_WAV_NORM_SPIN_LOWERED, this parameter
 *                           indicates which spin number the wavelets
 *                           were lowered from. Otherwise, it is ignored.
 * \retval none
 */
void s2let_wav_analysis_mw_real(
    double *f_wav,
    double *f_scal,
    const double *f,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;

    int verbosity = 0;
    ssht_dl_method_t dl_method = SSHT_DL_RISBO;
    so3_parameters_t so3_parameters = {};
    so3_parameters.verbosity = verbosity;
    so3_parameters.L = L;
    so3_parameters.N = N;
    so3_parameters.n_order = SO3_N_ORDER_NEGATIVE_FIRST;
    so3_parameters.storage = SO3_STORAGE_PADDED;
    so3_parameters.n_mode = SO3_N_MODE_ALL;
    so3_parameters.dl_method = dl_method;
    so3_parameters.reality = 1;

    int j, offset, offset_lmn;
    int J = s2let_j_max(parameters);
    //int l_min = s2let_axisym_el_min(B, J_min);

    complex double *wav_lm;
    double *scal_l;
    s2let_tiling_wavelet_allocate(&wav_lm, &scal_l, parameters);
    s2let_tiling_wavelet(wav_lm, scal_l, parameters);

    complex double *flm, *f_wav_lmn, *f_scal_lm;

    s2let_lm_allocate(&flm, L);
    ssht_core_mw_forward_sov_conv_sym_real(flm, f, L, dl_method, verbosity);

    s2let_allocate_f_wav_lmn(&f_wav_lmn, &f_scal_lm, parameters);
    s2let_wav_analysis_harmonic(f_wav_lmn, f_scal_lm, flm, wav_lm, scal_l, parameters);

    ssht_core_mw_inverse_sov_sym_real(f_scal, f_scal_lm, L, dl_method, verbosity);
    offset = 0;
    offset_lmn = 0;
    for (j = J_min; j <= J; ++j)
    {
        so3_core_inverse_via_ssht_real(
            f_wav + offset,
            f_wav_lmn + offset_lmn + (N-1) * L*L,
            &so3_parameters
        );
        offset_lmn += (2*N-1) * L*L;
        offset += (2*N-1) * L * (2*L-1);
    }

    free(flm);
    free(wav_lm);
    free(scal_l);
    free(f_scal_lm);
    free(f_wav_lmn);
}

/*!
 * Spherical wavelets : full resolution synthesis in real space, MW sampling.
 * Perform directional wavelet transform in real space (from scratch,
 * gives wavelet maps).
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Array of wavelets maps, MW sampling.
 * \param[in]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  normalization Indicates how to normalise the wavelets
 *                           and scaling function.
 * \param[in]  original_spin If normalization has value
 *                           S2LET_WAV_NORM_SPIN_LOWERED, this parameter
 *                           indicates which spin number the wavelets
 *                           were lowered from. Otherwise, it is ignored.
 * \retval none
 */
void s2let_wav_synthesis_mw_real(
    double *f,
    const double *f_wav,
    const double *f_scal,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;

    int verbosity = 0;
    ssht_dl_method_t dl_method = SSHT_DL_RISBO;
    so3_parameters_t so3_parameters = {};
    so3_parameters.verbosity = verbosity;
    so3_parameters.L = L;
    so3_parameters.N = N;
    so3_parameters.n_order = SO3_N_ORDER_NEGATIVE_FIRST;
    so3_parameters.storage = SO3_STORAGE_PADDED;
    so3_parameters.n_mode = SO3_N_MODE_ALL;
    so3_parameters.dl_method = dl_method;
    so3_parameters.reality = 1;

    int j, offset, offset_lmn;
    int J = s2let_j_max(parameters);
    //int l_min = s2let_axisym_el_min(B, J_min);

    complex double *wav_lm;
    double *scal_l;
    s2let_tiling_wavelet_allocate(&wav_lm, &scal_l, parameters);
    s2let_tiling_wavelet(wav_lm, scal_l, parameters);

    complex double *flm, *f_wav_lmn, *f_scal_lm;
    s2let_allocate_f_wav_lmn(&f_wav_lmn, &f_scal_lm, parameters);

    ssht_core_mw_forward_sov_conv_sym_real(f_scal_lm, f_scal, L, dl_method, verbosity);
    offset = 0;
    offset_lmn = 0;
    for (j = J_min; j <= J; ++j)
    {
        int n, el, m;

        so3_core_forward_via_ssht_real(
            f_wav_lmn + offset_lmn + (N-1) * L*L,
            f_wav + offset,
            &so3_parameters
        );

        // Expand lmn to negative n, because SO3 only
        // computes non-negative n
        for (n = 0; n < N; ++n)
        {
            for (el = 0; el < L; ++el)
            {
                for (m = -el; m <= el; ++m)
                {
                    int ind, nind;
                    so3_sampling_elmn2ind(&ind, el, m, n, &so3_parameters);
                    so3_sampling_elmn2ind(&nind, el, -m, -n, &so3_parameters);
                    int sign = (m+n)%2 ? -1 : 1;
                    f_wav_lmn[offset_lmn + nind] = sign * conj(f_wav_lmn[offset_lmn + ind]);
                }
            }
        }

        offset_lmn += (2*N-1) * L*L;
        offset += (2*N-1) * L * (2*L-1);
    }

    s2let_lm_allocate(&flm, L);
    s2let_wav_synthesis_harmonic(flm, f_wav_lmn, f_scal_lm, wav_lm, scal_l, parameters);

    ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

    free(flm);
    free(wav_lm);
    free(scal_l);
    free(f_scal_lm);
    free(f_wav_lmn);
}

/*!
 * Spherical wavelets : multi-resolution analysis in real space, MW sampling.
 * Perform directional wavelet transform in real space (from scratch,
 * gives wavelet maps).
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f_wav Array of wavelet maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  f Input function (MW sampling)
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  normalization Indicates how to normalise the wavelets
 *                           and scaling function.
 * \param[in]  original_spin If normalization has value
 *                           S2LET_WAV_NORM_SPIN_LOWERED, this parameter
 *                           indicates which spin number the wavelets
 *                           were lowered from. Otherwise, it is ignored.
 * \retval none
 */
void s2let_wav_analysis_mw_multires(
    complex double *f_wav,
    complex double *f_scal,
    const complex double *f,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;
    int spin = parameters->spin;

    int bandlimit;
    int verbosity = 0;
    ssht_dl_method_t dl_method = SSHT_DL_RISBO;
    so3_parameters_t so3_parameters = {};
    so3_parameters.verbosity = verbosity;
    so3_parameters.N = N;
    so3_parameters.n_order = SO3_N_ORDER_NEGATIVE_FIRST;
    so3_parameters.storage = SO3_STORAGE_PADDED;
    so3_parameters.n_mode = SO3_N_MODE_ALL;
    so3_parameters.dl_method = dl_method;

    int j, offset, offset_lmn;
    int J = s2let_j_max(parameters);
    //int l_min = s2let_axisym_el_min(B, J_min);

    complex double *wav_lm;
    double *scal_l;
    s2let_tiling_wavelet_allocate(&wav_lm, &scal_l, parameters);
    s2let_tiling_wavelet(wav_lm, scal_l, parameters);

    complex double *flm, *f_wav_lmn, *f_scal_lm;

    s2let_lm_allocate(&flm, L);
    ssht_core_mw_forward_sov_conv_sym(flm, f, L, spin, dl_method, verbosity);

    s2let_allocate_f_wav_lmn_multires(&f_wav_lmn, &f_scal_lm, parameters);
    s2let_wav_analysis_harmonic_multires(f_wav_lmn, f_scal_lm, flm, wav_lm, scal_l, parameters);

    bandlimit = MIN(s2let_bandlimit(J_min-1, parameters), L);
    // Note, this is a spin-0 transform!
    ssht_core_mw_inverse_sov_sym(f_scal, f_scal_lm, bandlimit, 0, dl_method, verbosity);
    offset = 0;
    offset_lmn = 0;
    for (j = J_min; j <= J; ++j)
    {
        bandlimit = MIN(s2let_bandlimit(j, parameters), L);
        so3_parameters.L = bandlimit;
        so3_core_inverse_via_ssht(
            f_wav + offset,
            f_wav_lmn + offset_lmn,
            &so3_parameters
        );
        offset_lmn += (2*N-1) * bandlimit*bandlimit;
        offset += (2*N-1) * bandlimit * (2*bandlimit-1);
    }

    free(flm);
    free(wav_lm);
    free(scal_l);
    free(f_scal_lm);
    free(f_wav_lmn);
}

/*!
 * Spherical wavelets : multi-resolution synthesis in real space, MW sampling.
 * Perform directional wavelet transform in real space (from scratch,
 * gives wavelet maps).
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Array of wavelets maps, MW sampling.
 * \param[in]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  normalization Indicates how to normalise the wavelets
 *                           and scaling function.
 * \param[in]  original_spin If normalization has value
 *                           S2LET_WAV_NORM_SPIN_LOWERED, this parameter
 *                           indicates which spin number the wavelets
 *                           were lowered from. Otherwise, it is ignored.
 * \retval none
 */
void s2let_wav_synthesis_mw_multires(
    complex double *f,
    const complex double *f_wav,
    const complex double *f_scal,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;
    int spin = parameters->spin;

    int bandlimit;
    int verbosity = 0;
    ssht_dl_method_t dl_method = SSHT_DL_RISBO;
    so3_parameters_t so3_parameters = {};
    so3_parameters.verbosity = verbosity;
    so3_parameters.N = N;
    so3_parameters.n_order = SO3_N_ORDER_NEGATIVE_FIRST;
    so3_parameters.storage = SO3_STORAGE_PADDED;
    so3_parameters.n_mode = SO3_N_MODE_ALL;
    so3_parameters.dl_method = dl_method;

    int j, offset, offset_lmn;
    int J = s2let_j_max(parameters);
    //int l_min = s2let_axisym_el_min(B, J_min);

    complex double *wav_lm;
    double *scal_l;
    s2let_tiling_wavelet_allocate(&wav_lm, &scal_l, parameters);
    s2let_tiling_wavelet(wav_lm, scal_l, parameters);

    complex double *flm, *f_wav_lmn, *f_scal_lm;
    s2let_allocate_f_wav_lmn_multires(&f_wav_lmn, &f_scal_lm, parameters);

    bandlimit = MIN(s2let_bandlimit(J_min-1, parameters), L);
    // Note, this is a spin-0 transform!
    ssht_core_mw_forward_sov_conv_sym(f_scal_lm, f_scal, bandlimit, 0, dl_method, verbosity);
    offset = 0;
    offset_lmn = 0;
    for (j = J_min; j <= J; ++j)
    {
        bandlimit = MIN(s2let_bandlimit(j, parameters), L);
        so3_parameters.L = bandlimit;
        so3_core_forward_via_ssht(
            f_wav_lmn + offset_lmn,
            f_wav + offset,
            &so3_parameters
        );

        offset_lmn += (2*N-1) * bandlimit*bandlimit;
        offset += (2*N-1) * bandlimit * (2*bandlimit-1);
    }

    s2let_lm_allocate(&flm, L);
    s2let_wav_synthesis_harmonic_multires(flm, f_wav_lmn, f_scal_lm, wav_lm, scal_l, parameters);

    ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

    free(flm);
    free(wav_lm);
    free(scal_l);
    free(f_scal_lm);
    free(f_wav_lmn);
}


/*!
 * Spherical wavelets : multi-resolution analysis in real space, MW sampling.
 * Perform directional wavelet transform in real space (from scratch,
 * gives wavelet maps).
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f_wav Array of wavelet maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  f Input function (MW sampling)
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  normalization Indicates how to normalise the wavelets
 *                           and scaling function.
 * \param[in]  original_spin If normalization has value
 *                           S2LET_WAV_NORM_SPIN_LOWERED, this parameter
 *                           indicates which spin number the wavelets
 *                           were lowered from. Otherwise, it is ignored.
 * \retval none
 */
void s2let_wav_analysis_mw_multires_real(
    double *f_wav,
    double *f_scal,
    const double *f,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;

    int bandlimit;
    int verbosity = 0;
    ssht_dl_method_t dl_method = SSHT_DL_RISBO;
    so3_parameters_t so3_parameters = {};
    so3_parameters.verbosity = verbosity;
    so3_parameters.N = N;
    so3_parameters.n_order = SO3_N_ORDER_NEGATIVE_FIRST;
    so3_parameters.storage = SO3_STORAGE_PADDED;
    so3_parameters.n_mode = SO3_N_MODE_ALL;
    so3_parameters.dl_method = dl_method;

    int j, offset, offset_lmn;
    int J = s2let_j_max(parameters);
    //int l_min = s2let_axisym_el_min(B, J_min);

    complex double *wav_lm;
    double *scal_l;
    s2let_tiling_wavelet_allocate(&wav_lm, &scal_l, parameters);
    s2let_tiling_wavelet(wav_lm, scal_l, parameters);

    complex double *flm, *f_wav_lmn, *f_scal_lm;

    s2let_lm_allocate(&flm, L);
    ssht_core_mw_forward_sov_conv_sym_real(flm, f, L, dl_method, verbosity);

    s2let_allocate_f_wav_lmn_multires(&f_wav_lmn, &f_scal_lm, parameters);
    s2let_wav_analysis_harmonic_multires(f_wav_lmn, f_scal_lm, flm, wav_lm, scal_l, parameters);

    bandlimit = MIN(s2let_bandlimit(J_min-1, parameters), L);
    // Note, this is a spin-0 transform!
    ssht_core_mw_inverse_sov_sym_real(f_scal, f_scal_lm, bandlimit, dl_method, verbosity);
    offset = 0;
    offset_lmn = 0;
    for (j = J_min; j <= J; ++j)
    {
        bandlimit = MIN(s2let_bandlimit(j, parameters), L);
        so3_parameters.L = bandlimit;
        so3_core_inverse_via_ssht_real(
            f_wav + offset,
            f_wav_lmn + offset_lmn + (N-1) * bandlimit*bandlimit,
            &so3_parameters
        );
        offset_lmn += (2*N-1) * bandlimit*bandlimit;
        offset += (2*N-1) * bandlimit * (2*bandlimit-1);
    }

    free(flm);
    free(wav_lm);
    free(scal_l);
    free(f_scal_lm);
    free(f_wav_lmn);
}

/*!
 * Spherical wavelets : multi-resolution synthesis in real space, MW sampling.
 * Perform directional wavelet transform in real space (from scratch,
 * gives wavelet maps).
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Array of wavelets maps, MW sampling.
 * \param[in]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  normalization Indicates how to normalise the wavelets
 *                           and scaling function.
 * \param[in]  original_spin If normalization has value
 *                           S2LET_WAV_NORM_SPIN_LOWERED, this parameter
 *                           indicates which spin number the wavelets
 *                           were lowered from. Otherwise, it is ignored.
 * \retval none
 */
void s2let_wav_synthesis_mw_multires_real(
    double *f,
    const double *f_wav,
    const double *f_scal,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;

    int bandlimit;
    int verbosity = 0;
    ssht_dl_method_t dl_method = SSHT_DL_RISBO;
    so3_parameters_t so3_parameters = {};
    so3_parameters.verbosity = verbosity;
    so3_parameters.N = N;
    so3_parameters.n_order = SO3_N_ORDER_NEGATIVE_FIRST;
    so3_parameters.storage = SO3_STORAGE_PADDED;
    so3_parameters.n_mode = SO3_N_MODE_ALL;
    so3_parameters.dl_method = dl_method;

    int j, offset, offset_lmn;
    int J = s2let_j_max(parameters);
    //int l_min = s2let_axisym_el_min(B, J_min);

    complex double *wav_lm;
    double *scal_l;
    s2let_tiling_wavelet_allocate(&wav_lm, &scal_l, parameters);
    s2let_tiling_wavelet(wav_lm, scal_l, parameters);

    complex double *flm, *f_wav_lmn, *f_scal_lm;
    s2let_allocate_f_wav_lmn_multires(&f_wav_lmn, &f_scal_lm, parameters);

    bandlimit = MIN(s2let_bandlimit(J_min-1, parameters), L);
    // Note, this is a spin-0 transform!
    ssht_core_mw_forward_sov_conv_sym_real(f_scal_lm, f_scal, bandlimit, dl_method, verbosity);
    offset = 0;
    offset_lmn = 0;
    for (j = J_min; j <= J; ++j)
    {
        int n, el, m;

        bandlimit = MIN(s2let_bandlimit(j, parameters), L);
        so3_parameters.L = bandlimit;
        so3_core_forward_via_ssht_real(
            f_wav_lmn + offset_lmn + (N-1) * bandlimit*bandlimit,
            f_wav + offset,
            &so3_parameters
        );

        // Expand lmn to negative n, because SO3 only
        // computes non-negative n
        for (n = 0; n < N; ++n)
        {
            for (el = 0; el < bandlimit; ++el)
            {
                for (m = -el; m <= el; ++m)
                {
                    int ind, nind;
                    so3_sampling_elmn2ind(&ind, el, m, n, &so3_parameters);
                    so3_sampling_elmn2ind(&nind, el, -m, -n, &so3_parameters);
                    int sign = (m+n)%2 ? -1 : 1;
                    f_wav_lmn[offset_lmn + nind] = sign * conj(f_wav_lmn[offset_lmn + ind]);
                }
            }
        }

        offset_lmn += (2*N-1) * bandlimit*bandlimit;
        offset += (2*N-1) * bandlimit * (2*bandlimit-1);
    }

    s2let_lm_allocate(&flm, L);
    s2let_wav_synthesis_harmonic_multires(flm, f_wav_lmn, f_scal_lm, wav_lm, scal_l, parameters);

    ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

    free(flm);
    free(wav_lm);
    free(scal_l);
    free(f_scal_lm);
    free(f_wav_lmn);
}

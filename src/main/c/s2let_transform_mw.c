// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <complex.h>
#include <ssht.h>
#include <so3.h>
#include <stdlib.h>

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
    int spin = parameters->spin;
    ssht_dl_method_t dl_method = parameters->dl_method;
    int verbosity = parameters->verbosity;

    complex double *flm;
    s2let_lm_allocate(&flm, L);

    switch (parameters->sampling_scheme)
    {
    case S2LET_SAMPLING_MW:
        ssht_core_mw_forward_sov_conv_sym(flm, f, L, spin, dl_method, verbosity);
        break;
    case S2LET_SAMPLING_MW_SS:
        ssht_core_mw_forward_sov_conv_sym_ss(flm, f, L, spin, dl_method, verbosity);
        break;
    default:
        S2LET_ERROR_GENERIC("Sampling scheme not supported.");
    }

    s2let_wav_analysis_lm2wav(f_wav, f_scal, flm, parameters);

    free(flm);
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
    int spin = parameters->spin;
    ssht_dl_method_t dl_method = parameters->dl_method;
    int verbosity = 0;

    complex double *flm;
    s2let_lm_allocate(&flm, L);

    s2let_wav_synthesis_lm2wav(flm, f_wav, f_scal, parameters);

    switch (parameters->sampling_scheme)
    {
    case S2LET_SAMPLING_MW:
        ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);
        break;
    case S2LET_SAMPLING_MW_SS:
        ssht_core_mw_inverse_sov_sym_ss(f, flm, L, spin, dl_method, verbosity);
        break;
    default:
        S2LET_ERROR_GENERIC("Sampling scheme not supported.");
    }

    free(flm);
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
    ssht_dl_method_t dl_method = parameters->dl_method;
    int verbosity = 0;

    complex double *flm;
    s2let_lm_allocate(&flm, L);

    switch (parameters->sampling_scheme)
    {
    case S2LET_SAMPLING_MW:
        ssht_core_mw_forward_sov_conv_sym_real(flm, f, L, dl_method, verbosity);
        break;
    case S2LET_SAMPLING_MW_SS:
        ssht_core_mw_forward_sov_conv_sym_ss_real(flm, f, L, dl_method, verbosity);
        break;
    default:
        S2LET_ERROR_GENERIC("Sampling scheme not supported.");
    }

    s2let_wav_analysis_lm2wav_real(f_wav, f_scal, flm, parameters);

    free(flm);
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
    ssht_dl_method_t dl_method = parameters->dl_method;
    int verbosity = 0;

    complex double *flm;
    s2let_lm_allocate(&flm, L);

    s2let_wav_synthesis_lm2wav_real(flm, f_wav, f_scal, parameters);

    switch (parameters->sampling_scheme)
    {
    case S2LET_SAMPLING_MW:
        ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);
        break;
    case S2LET_SAMPLING_MW_SS:
        ssht_core_mw_inverse_sov_sym_ss_real(f, flm, L, dl_method, verbosity);
        break;
    default:
        S2LET_ERROR_GENERIC("Sampling scheme not supported.");
    }

    free(flm);
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
    int spin = parameters->spin;
    ssht_dl_method_t dl_method = parameters->dl_method;
    int verbosity = 0;

    complex double *flm;
    s2let_lm_allocate(&flm, L);

    switch (parameters->sampling_scheme)
    {
    case S2LET_SAMPLING_MW:
        ssht_core_mw_forward_sov_conv_sym(flm, f, L, spin, dl_method, verbosity);
        break;
    case S2LET_SAMPLING_MW_SS:
        ssht_core_mw_forward_sov_conv_sym_ss(flm, f, L, spin, dl_method, verbosity);
        break;
    default:
        S2LET_ERROR_GENERIC("Sampling scheme not supported.");
    }

    s2let_wav_analysis_lm2wav_multires(f_wav, f_scal, flm, parameters);

    free(flm);
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
    int spin = parameters->spin;
    ssht_dl_method_t dl_method = parameters->dl_method;
    int verbosity = 0;

    complex double *flm;
    s2let_lm_allocate(&flm, L);

    s2let_wav_synthesis_lm2wav_multires(flm, f_wav, f_scal, parameters);

    switch (parameters->sampling_scheme)
    {
    case S2LET_SAMPLING_MW:
        ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);
        break;
    case S2LET_SAMPLING_MW_SS:
        ssht_core_mw_inverse_sov_sym_ss(f, flm, L, spin, dl_method, verbosity);
        break;
    default:
        S2LET_ERROR_GENERIC("Sampling scheme not supported.");
    }

    free(flm);
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
    ssht_dl_method_t dl_method = parameters->dl_method;
    int verbosity = 0;

    complex double *flm;
    s2let_lm_allocate(&flm, L);

    switch (parameters->sampling_scheme)
    {
    case S2LET_SAMPLING_MW:
        ssht_core_mw_forward_sov_conv_sym_real(flm, f, L, dl_method, verbosity);
        break;
    case S2LET_SAMPLING_MW_SS:
        ssht_core_mw_forward_sov_conv_sym_ss_real(flm, f, L, dl_method, verbosity);
        break;
    default:
        S2LET_ERROR_GENERIC("Sampling scheme not supported.");
    }

    s2let_wav_analysis_lm2wav_multires_real(f_wav, f_scal, flm, parameters);

    free(flm);
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
    ssht_dl_method_t dl_method = parameters->dl_method;
    int verbosity = 0;

    complex double *flm;
    s2let_lm_allocate(&flm, L);

    s2let_wav_synthesis_lm2wav_multires_real(flm, f_wav, f_scal, parameters);

    switch (parameters->sampling_scheme)
    {
    case S2LET_SAMPLING_MW:
        ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);
        break;
    case S2LET_SAMPLING_MW_SS:
        ssht_core_mw_inverse_sov_sym_ss_real(f, flm, L, dl_method, verbosity);
        break;
    default:
        S2LET_ERROR_GENERIC("Sampling scheme not supported.");
    }

    free(flm);
}

// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <ssht.h>
#include <assert.h>
#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*!
 * Test the identity relation of the wavelet tiling in harmonic space.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_tiling_axisym_test(int B, int L, int J_min)
{
  double *kappa, *kappa0;

  // Allocate the kernels corresponding to the parameters B, L
  s2let_tiling_axisym_allocate(&kappa, &kappa0, B, L);

  // Construct the tiling of harmonic space
  s2let_tiling_axisym(kappa, kappa0, B, L, J_min);

  // Check that they recover the identity relation,
  // ensuring exactness of the wavelet transform.
  double res = s2let_tiling_axisym_check_identity(kappa, kappa0, B, L, J_min);
  printf("  - Maximum error : %6.5e\n", res);

  free(kappa);
  free(kappa0);
}

/*!
 * Test the identity relation of the directionality components for
 * the wavelet tiling in harmonic space.
 *
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Azimuthal band-limit.
 * \retval none
 */
void s2let_tiling_direction_test(int L, int N)
{
  complex double *s_elm;
  double error;

  // Allocate space for the harmonic coefficients
  s2let_tiling_direction_allocate(&s_elm, L, N);

  // Construct the harmonic coefficients
  s2let_tiling_direction(s_elm, L, N);

  // Check that they recover the identity relation,
  // ensuring exactness of the wavelet transform.
  error = s2let_tiling_direction_check_identity(s_elm, L, N);
  printf("  - Maximum error : %6.5e\n", error);

  free(s_elm);
}

/*!
 * Test the identity relation of the directional wavelets for
 * the wavelet tiling in harmonic space.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \retval none
 */
void s2let_tiling_wavelet_test(int B, int L, int J_min, int N, int spin)
{
  complex double *phi;
  double *psi;
  double error;

  // Allocate space for the harmonic coefficients
  s2let_tiling_wavelet_allocate(&phi, &psi, B, L, N);

  // Construct the harmonic coefficients
  s2let_tiling_wavelet(phi, psi, B, L, J_min, N, spin, S2LET_WAV_NORM_DEFAULT, 0);

  // Check that they recover the identity relation,
  // ensuring exactness of the wavelet transform.
  error = s2let_tiling_wavelet_check_identity(phi, psi, B, L, J_min, N, spin);
  printf("  - Maximum error : %6.5e\n", error);

  free(phi);
  free(psi);
}

void s2let_binomial_coefficient_test(int n_max)
{
    const int nRepeat = 100000;
    int n, k, i;
    int firstError = 0;

    long error;

    clock_t time_start, time_end;

    for (n = 1; n <= n_max; ++n)
    {
        for (k = 0; k <= n/2; ++k)
        {
            error = binomial_coefficient(n, k, 0) -
                        binomial_coefficient(n, k, 1);

            if (error && !firstError)
            {
                printf("  - First error at: n = %d, k = %d, error = %ld\n", n, k, error);
                firstError = 1;
            }
        }
    }

    printf("  - Maximum error: %ld\n", error);

    printf("  - Duration for %d computations of (%d,%d)\n", nRepeat, n_max, n_max/2);

    time_start = clock();
    for (i = 0; i < nRepeat; ++i)
        binomial_coefficient(n_max, n_max/2, 0);
    time_end = clock();
    printf("    logfact implementation: %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    time_start = clock();
    for (i = 0; i < nRepeat; ++i)
        binomial_coefficient(n_max, n_max/2, 1);
    time_end = clock();
    printf("    exact implementation: %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
}

/*!
 * Test the exactness of the full resolution wavelet transform in harmonic space.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_transform_axisym_lm_wav_test(int B, int L, int J_min, int seed)
{
  clock_t time_start, time_end;
  double *wav_lm, *scal_lm;

  // Allocate the wavelet kernels
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, B, L);

  // Compute the wavelet kernels
  time_start = clock();
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, B, L, J_min);
  time_end = clock();
  printf("  - Generate wavelets  : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  complex double *f_wav_lm, *f_scal_lm, *flm, *flm_rec;
  s2let_lm_allocate(&flm, L);
  s2let_lm_allocate(&flm_rec, L);

  // Generate a random spherical harmonic decomposition
  s2let_lm_random_flm(flm, L, 0, seed);

  // Allocate space for the wavelet scales (their harmonic coefficients)
  s2let_transform_axisym_lm_allocate_f_wav(&f_wav_lm, &f_scal_lm, B, L, J_min);

  // Perform the wavelet transform through exact harmonic tiling
  time_start = clock();
  s2let_transform_axisym_lm_wav_analysis(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, B, L, J_min);
  time_end = clock();
  printf("  - Wavelet analysis   : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Reconstruct the initial harmonic coefficients from those of the wavelets
  time_start = clock();
  s2let_transform_axisym_lm_wav_synthesis(flm_rec, f_wav_lm, f_scal_lm, wav_lm, scal_lm, B, L, J_min);
  time_end = clock();
  printf("  - Wavelet synthesis  : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Compute the maximum absolute error on the harmonic coefficients
  printf("  - Maximum abs error  : %6.5e\n",
	 maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

  free(flm);
  free(flm_rec);
  free(f_wav_lm);
  free(f_scal_lm);
  free(wav_lm);
  free(scal_lm);
}

/*!
 * Test the exactness of the multiresolution wavelet transform in harmonic space.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_transform_axisym_lm_wav_multires_test(int B, int L, int J_min, int seed)
{
  clock_t time_start, time_end;
  double *wav_lm, *scal_lm;

  // Allocate the wavelet kernels
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, B, L);

  // Compute the wavelet kernels
  time_start = clock();
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, B, L, J_min);
  time_end = clock();
  printf("  - Generate wavelets  : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  complex double *f_wav_lm, *f_scal_lm, *flm, *flm_rec;
  s2let_lm_allocate(&flm, L);
  s2let_lm_allocate(&flm_rec, L);

  // Generate a random spherical harmonic decomposition
  s2let_lm_random_flm(flm, L, 0, seed);

  // Allocate space for the wavelet scales (their harmonic coefficients)
  s2let_transform_axisym_lm_allocate_f_wav_multires(&f_wav_lm, &f_scal_lm, B, L, J_min);

  // Perform the wavelet transform through exact harmonic tiling
  time_start = clock();
  s2let_transform_axisym_lm_wav_analysis_multires(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, B, L, J_min);
  time_end = clock();
  printf("  - Wavelet analysis   : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Reconstruct the initial harmonic coefficients from those of the wavelets
  time_start = clock();
  s2let_transform_axisym_lm_wav_synthesis_multires(flm_rec, f_wav_lm, f_scal_lm, wav_lm, scal_lm, B, L, J_min);
  time_end = clock();
  printf("  - Wavelet synthesis  : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Compute the maximum absolute error on the harmonic coefficients
  printf("  - Maximum abs error  : %6.5e\n",
	 maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

  free(flm);
  free(flm_rec);
  free(f_wav_lm);
  free(f_scal_lm);
  free(wav_lm);
  free(scal_lm);
}

/*!
 * Test the exactness of the full resolution directional wavelet transform
 * in harmonic space.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_harmonic_test(int B, int L, int J_min, int N, int spin, int seed)
{
  clock_t time_start, time_end;
  complex double *psi;
  double *phi;

  // Allocate the wavelet kernels
  s2let_tiling_wavelet_allocate(&psi, &phi, B, L, N);

  // Compute the wavelet kernels
  time_start = clock();
  s2let_tiling_wavelet(psi, phi, B, L, J_min, N, spin, S2LET_WAV_NORM_DEFAULT, 0);
  time_end = clock();
  printf("  - Generate wavelets  : %4.4f seconds\n",
     (time_end - time_start) / (double)CLOCKS_PER_SEC);

  complex double *f_wav_lmn, *f_scal_lm, *flm, *flm_rec;
  s2let_lm_allocate(&flm, L);
  s2let_lm_allocate(&flm_rec, L);

  // Generate a random spherical harmonic decomposition
  s2let_lm_random_flm(flm, L, spin, seed);

  // Allocate space for the wavelet scales (their harmonic/Wigner coefficients)
  s2let_allocate_f_wav_lmn(&f_wav_lmn, &f_scal_lm, B, L, J_min, N);

  // Perform the wavelet transform through exact harmonic tiling
  time_start = clock();
  s2let_wav_analysis_harmonic(f_wav_lmn, f_scal_lm, flm, psi, phi, B, L, J_min, N, spin);
  time_end = clock();
  printf("  - Wavelet analysis   : %4.4f seconds\n",
     (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Reconstruct the initial harmonic coefficients from those of the wavelets
  time_start = clock();
  s2let_wav_synthesis_harmonic(flm_rec, f_wav_lmn, f_scal_lm, psi, phi, B, L, J_min, N, spin);
  time_end = clock();
  printf("  - Wavelet synthesis  : %4.4f seconds\n",
     (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Compute the maximum absolute error on the harmonic coefficients
  printf("  - Maximum abs error  : %6.5e\n",
     maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

  free(flm);
  free(flm_rec);
  free(f_wav_lmn);
  free(f_scal_lm);
  free(psi);
  free(phi);
}

/*!
 * Test the exactness of the multi-resolution directional wavelet transform
 * in harmonic space.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_harmonic_multires_test(int B, int L, int J_min, int N, int spin, int seed)
{
  clock_t time_start, time_end;
  complex double *psi;
  double *phi;

  // Allocate the wavelet kernels
  s2let_tiling_wavelet_allocate(&psi, &phi, B, L, N);

  // Compute the wavelet kernels
  time_start = clock();
  s2let_tiling_wavelet(psi, phi, B, L, J_min, N, spin, S2LET_WAV_NORM_DEFAULT, 0);
  time_end = clock();
  printf("  - Generate wavelets  : %4.4f seconds\n",
     (time_end - time_start) / (double)CLOCKS_PER_SEC);

  complex double *f_wav_lmn, *f_scal_lm, *flm, *flm_rec;
  s2let_lm_allocate(&flm, L);
  s2let_lm_allocate(&flm_rec, L);

  // Generate a random spherical harmonic decomposition
  s2let_lm_random_flm(flm, L, spin, seed);

  // Allocate space for the wavelet scales (their harmonic/Wigner coefficients)
  s2let_allocate_f_wav_lmn_multires(&f_wav_lmn, &f_scal_lm, B, L, J_min, N);

  // Perform the wavelet transform through exact harmonic tiling
  time_start = clock();
  s2let_wav_analysis_harmonic_multires(f_wav_lmn, f_scal_lm, flm, psi, phi, B, L, J_min, N, spin);
  time_end = clock();
  printf("  - Wavelet analysis   : %4.4f seconds\n",
     (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Reconstruct the initial harmonic coefficients from those of the wavelets
  time_start = clock();
  s2let_wav_synthesis_harmonic_multires(flm_rec, f_wav_lmn, f_scal_lm, psi, phi, B, L, J_min, N, spin);
  time_end = clock();
  printf("  - Wavelet synthesis  : %4.4f seconds\n",
     (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Compute the maximum absolute error on the harmonic coefficients
  printf("  - Maximum abs error  : %6.5e\n",
     maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

  free(flm);
  free(flm_rec);
  free(f_wav_lmn);
  free(f_scal_lm);
  free(psi);
  free(phi);
}


/*!
 * Test the exactness of the full resolution wavelet transform in real space for complex functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_transform_axisym_wav_test(int B, int L, int J_min, int seed)
{
  clock_t time_start, time_end;
  int spin = 0;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;
  //int J = s2let_j_max(L, B);

  complex double *f, *f_rec, *flm, *flm_rec;
  s2let_lm_allocate(&flm, L);
  s2let_lm_allocate(&flm_rec, L);
  s2let_mw_allocate(&f, L);
  s2let_mw_allocate(&f_rec, L);

  // Generate random harmonic coefficients for a complex signal
  s2let_lm_random_flm(flm, L, 0, seed);

  // Construct the corresponding signal on the sphere (MW sampling)
  ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

  // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
  complex double *f_wav, *f_scal;
  s2let_transform_axisym_allocate_mw_f_wav(&f_wav, &f_scal, B, L, J_min);

  // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
  time_start = clock();
  s2let_transform_axisym_wav_analysis_mw(f_wav, f_scal, f, B, L, J_min);
  time_end = clock();
  printf("  - Wavelet analysis   : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Reconstruct the initial signal from the wavelet maps from scratch
  time_start = clock();
  s2let_transform_axisym_wav_synthesis_mw(f_rec, f_wav, f_scal, B, L, J_min);
  time_end = clock();
  printf("  - Wavelet synthesis  : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Compute the initial harmonic coefficients back
  ssht_core_mw_forward_sov_conv_sym(flm_rec, f_rec, L, spin, dl_method, verbosity);

  // Compute the maximum absolute error on the harmonic coefficients
  printf("  - Maximum abs error  : %6.5e\n",
	 maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

  free(f);
  free(f_rec);
  free(f_wav);
  free(f_scal);
}

/*!
 * Test the exactness of the full resolution wavelet transform in real space for real functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_transform_axisym_wav_real_test(int B, int L, int J_min, int seed)
{
  clock_t time_start, time_end;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;
  //int J = s2let_j_max(L, B);

  complex *flm, *flm_rec;
  double *f, *f_rec;
  s2let_lm_allocate(&flm, L);
  s2let_lm_allocate(&flm_rec, L);
  s2let_mw_allocate_real(&f, L);
  s2let_mw_allocate_real(&f_rec, L);

  // Generate random harmonic coefficients for a real signal
  s2let_lm_random_flm_real(flm, L, seed);

  // Construct the corresponding signal on the sphere (MW sampling)
  ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

  // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
  double *f_wav, *f_scal;
  s2let_transform_axisym_allocate_mw_f_wav_real(&f_wav, &f_scal, B, L, J_min);

  // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
  time_start = clock();
  s2let_transform_axisym_wav_analysis_mw_real(f_wav, f_scal, f, B, L, J_min);
  time_end = clock();
  printf("  - Wavelet analysis   : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Reconstruct the initial signal from the wavelet maps from scratch
  time_start = clock();
  s2let_transform_axisym_wav_synthesis_mw_real(f_rec, f_wav, f_scal, B, L, J_min);
  time_end = clock();
  printf("  - Wavelet synthesis  : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Compute the initial harmonic coefficients back
  ssht_core_mw_forward_sov_conv_sym_real(flm_rec, f_rec, L, dl_method, verbosity);

  // Compute the maximum absolute error on the harmonic coefficients
  printf("  - Maximum abs error  : %6.5e\n",
	 maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

  free(f);
  free(f_rec);
  free(f_wav);
  free(f_scal);
}


/*!
 * Test the exactness of the multiresolution wavelet transform in real space for complex functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_transform_axisym_wav_multires_test(int B, int L, int J_min, int seed)
{
  clock_t time_start, time_end;
  int spin = 0;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;
  //int J = s2let_j_max(L, B);

  complex double *f, *f_rec, *flm, *flm_rec;
  s2let_lm_allocate(&flm, L);
  s2let_lm_allocate(&flm_rec, L);
  s2let_mw_allocate(&f, L);
  s2let_mw_allocate(&f_rec, L);

  // Generate random harmonic coefficients for a complex signal
  s2let_lm_random_flm(flm, L, 0, seed);

  // Construct the corresponding signal on the sphere (MW sampling)
  ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

  // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
  complex double *f_wav, *f_scal;
  s2let_transform_axisym_allocate_mw_f_wav_multires(&f_wav, &f_scal, B, L, J_min);

  // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
  time_start = clock();
  s2let_transform_axisym_wav_analysis_mw_multires(f_wav, f_scal, f, B, L, J_min);
  time_end = clock();
  printf("  - Wavelet analysis   : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Reconstruct the initial signal from the wavelet maps from scratch
  time_start = clock();
  s2let_transform_axisym_wav_synthesis_mw_multires(f_rec, f_wav, f_scal, B, L, J_min);
  time_end = clock();
  printf("  - Wavelet synthesis  : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Compute the initial harmonic coefficients back
  ssht_core_mw_forward_sov_conv_sym(flm_rec, f_rec, L, spin, dl_method, verbosity);

  // Compute the maximum absolute error on the harmonic coefficients
  printf("  - Maximum abs error  : %6.5e\n",
	 maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

  free(f);
  free(f_rec);
  free(f_wav);
  free(f_scal);
}


/*!
 * Test the exactness of the multiresolution wavelet transform in real space for real functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_transform_axisym_wav_multires_real_test(int B, int L, int J_min, int seed)
{
  clock_t time_start, time_end;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;
  //int J = s2let_j_max(L, B);

  complex *flm, *flm_rec;
  double *f, *f_rec;
  s2let_lm_allocate(&flm, L);
  s2let_lm_allocate(&flm_rec, L);
  s2let_mw_allocate_real(&f, L);
  s2let_mw_allocate_real(&f_rec, L);

  // Generate random harmonic coefficients for a real signal
  s2let_lm_random_flm_real(flm, L, seed);

  // Construct the corresponding signal on the sphere (MW sampling)
  ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

  // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
  double *f_wav, *f_scal;
  s2let_transform_axisym_allocate_mw_f_wav_multires_real(&f_wav, &f_scal, B, L, J_min);

  // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
  time_start = clock();
  s2let_transform_axisym_wav_analysis_mw_multires_real(f_wav, f_scal, f, B, L, J_min);
  time_end = clock();
  printf("  - Wavelet analysis   : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Reconstruct the initial signal from the wavelet maps from scratch
  time_start = clock();
  s2let_transform_axisym_wav_synthesis_mw_multires_real(f_rec, f_wav, f_scal, B, L, J_min);
  time_end = clock();
  printf("  - Wavelet synthesis  : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Compute the initial harmonic coefficients back
  ssht_core_mw_forward_sov_conv_sym_real(flm_rec, f_rec, L, dl_method, verbosity);

  // Compute the maximum absolute error on the harmonic coefficients
  printf("  - Maximum abs error  : %6.5e\n",
	 maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

  free(f);
  free(f_rec);
  free(f_wav);
  free(f_scal);
}



/*!
 * Test the exactness of the full resolution directional wavelet transform
 * in pixel space for complex functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_mw_test(int B, int L, int J_min, int N, int spin, int seed)
{
    clock_t time_start, time_end;

    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = spin;

    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);

    complex double *f, *f_rec, *flm, *flm_rec;
    s2let_lm_allocate(&flm, L);
    s2let_lm_allocate(&flm_rec, L);
    s2let_mw_allocate(&f, L);
    s2let_mw_allocate(&f_rec, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);

    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    complex double *f_wav, *f_scal;
    s2let_allocate_mw_f_wav(&f_wav, &f_scal, &parameters);

    // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_wav_analysis_mw(f_wav, f_scal, f, B, L, J_min, N, spin, S2LET_WAV_NORM_DEFAULT, 0);
    time_end = clock();
    printf("  - Wavelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Reconstruct the initial signal from the wavelet maps from scratch
    time_start = clock();
    s2let_wav_synthesis_mw(f_rec, f_wav, f_scal, B, L, J_min, N, spin, S2LET_WAV_NORM_DEFAULT, 0);
    time_end = clock();
    printf("  - Wavelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym(flm_rec, f_rec, L, spin, dl_method, verbosity);

    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_wav);
    free(f_scal);
}

/*!
 * Test the exactness of the full resolution directional wavelet transform
 * in pixel space for real functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_mw_real_test(int B, int L, int J_min, int N, int seed)
{
    clock_t time_start, time_end;

    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = 0;

    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);

    double *f, *f_rec;
    complex double *flm, *flm_rec;
    s2let_lm_allocate(&flm, L);
    s2let_lm_allocate(&flm_rec, L);
    s2let_mw_allocate_real(&f, L);
    s2let_mw_allocate_real(&f_rec, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm_real(flm, L, seed);

    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    double *f_wav, *f_scal;
    s2let_allocate_mw_f_wav_real(&f_wav, &f_scal, &parameters);

    // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_wav_analysis_mw_real(f_wav, f_scal, f, B, L, J_min, N);
    time_end = clock();
    printf("  - Wavelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Reconstruct the initial signal from the wavelet maps from scratch
    time_start = clock();
    s2let_wav_synthesis_mw_real(f_rec, f_wav, f_scal, B, L, J_min, N);
    time_end = clock();
    printf("  - Wavelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym_real(flm_rec, f_rec, L, dl_method, verbosity);

    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_wav);
    free(f_scal);
}

/*!
 * Test the exactness of the multi-resolution directional wavelet transform
 * in pixel space for complex functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_mw_multires_test(int B, int L, int J_min, int N, int spin, int seed)
{
    clock_t time_start, time_end;

    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = spin;

    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);

    complex double *f, *f_rec, *flm, *flm_rec;
    s2let_lm_allocate(&flm, L);
    s2let_lm_allocate(&flm_rec, L);
    s2let_mw_allocate(&f, L);
    s2let_mw_allocate(&f_rec, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);

    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    complex double *f_wav, *f_scal;
    s2let_allocate_mw_f_wav_multires(&f_wav, &f_scal, &parameters);

    // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_wav_analysis_mw_multires(f_wav, f_scal, f, B, L, J_min, N, spin, S2LET_WAV_NORM_DEFAULT, 0);
    time_end = clock();
    printf("  - Wavelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Reconstruct the initial signal from the wavelet maps from scratch
    time_start = clock();
    s2let_wav_synthesis_mw_multires(f_rec, f_wav, f_scal, B, L, J_min, N, spin, S2LET_WAV_NORM_DEFAULT, 0);
    time_end = clock();
    printf("  - Wavelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym(flm_rec, f_rec, L, spin, dl_method, verbosity);

    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_wav);
    free(f_scal);
}

/*!
 * Test the exactness of the multi-resolution directional wavelet transform
 * in pixel space for real functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_mw_multires_real_test(int B, int L, int J_min, int N, int seed)
{
    clock_t time_start, time_end;

    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = 0;

    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);

    double *f, *f_rec;
    complex double *flm, *flm_rec;
    s2let_lm_allocate(&flm, L);
    s2let_lm_allocate(&flm_rec, L);
    s2let_mw_allocate_real(&f, L);
    s2let_mw_allocate_real(&f_rec, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm_real(flm, L, seed);

    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    double *f_wav, *f_scal;
    s2let_allocate_mw_f_wav_multires_real(&f_wav, &f_scal, &parameters);

    // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_wav_analysis_mw_multires_real(f_wav, f_scal, f, B, L, J_min, N);
    time_end = clock();
    printf("  - Wavelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Reconstruct the initial signal from the wavelet maps from scratch
    time_start = clock();
    s2let_wav_synthesis_mw_multires_real(f_rec, f_wav, f_scal, B, L, J_min, N);
    time_end = clock();
    printf("  - Wavelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym_real(flm_rec, f_rec, L, dl_method, verbosity);

    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_wav);
    free(f_scal);
}

/*!
 * Test that the directional algorithms reduce to the axisymmetric ones,
 * provided spin = 0, N = 1
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_transform_axisym_vs_directional_mw_test(B, L, J_min, seed)
{
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;

    int spin = parameters.spin = 0;
    int N = parameters.N = 1;
    int J = s2let_j_max(&parameters);
    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;

    int i;

    double wav_error, scal_error;

    complex double *f, *flm;
    s2let_lm_allocate(&flm, L);
    s2let_mw_allocate(&f, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);

    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    // from both transforms.
    complex double *f_wav_axisym, *f_scal_axisym, *f_wav_dir, *f_scal_dir;
    s2let_transform_axisym_allocate_mw_f_wav(&f_wav_axisym, &f_scal_axisym, B, L, J_min);
    s2let_allocate_mw_f_wav(&f_wav_dir, &f_scal_dir, &parameters);

    // Do both transforms
    s2let_transform_axisym_wav_analysis_mw(f_wav_axisym, f_scal_axisym, f, B, L, J_min);
    s2let_wav_analysis_mw(f_wav_dir, f_scal_dir, f, B, L, J_min, N, spin, S2LET_WAV_NORM_DEFAULT, 0);

    // Account for the different wavelet normalisations
    for (i = 0; i < (J-J_min+1)*L*(2*L-1); ++i)
        f_wav_dir[i] *= sqrt(2*PI);

    // Compute the maximum absolute error in the computed wavelet transform
    wav_error = maxerr_cplx(f_wav_axisym, f_wav_dir, (J-J_min+1)*L*(2*L-1));
    scal_error = maxerr_cplx(f_scal_axisym, f_scal_dir, L*(2*L-1));

    printf("  - Maximum abs error in wavelets :         %6.5e\n", wav_error);
    printf("  - Maximum abs error in scaling function : %6.5e\n", scal_error);
    fflush(NULL);
}

/*!
 * Test that the directional multi-resolution algorithms reduce to
 * the axisymmetric ones, provided spin = 0, N = 1
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_transform_axisym_vs_directional_mw_multires_test(B, L, J_min, seed)
{
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;

    int spin = parameters.spin = 0;
    int N = parameters.N = 1;
    int J = s2let_j_max(&parameters);
    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;

    int samples, bandlimit, i, j;

    double wav_error, scal_error;

    complex double *f, *flm;
    s2let_lm_allocate(&flm, L);
    s2let_mw_allocate(&f, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);

    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    // from both transforms.
    complex double *f_wav_axisym, *f_scal_axisym, *f_wav_dir, *f_scal_dir;
    s2let_transform_axisym_allocate_mw_f_wav_multires(&f_wav_axisym, &f_scal_axisym, B, L, J_min);
    s2let_allocate_mw_f_wav_multires(&f_wav_dir, &f_scal_dir, &parameters);

    // Do both transforms
    s2let_transform_axisym_wav_analysis_mw_multires(f_wav_axisym, f_scal_axisym, f, B, L, J_min);
    s2let_wav_analysis_mw_multires(f_wav_dir, f_scal_dir, f, B, L, J_min, N, spin, S2LET_WAV_NORM_DEFAULT, 0);

    samples = 0;
    for (j = J_min; j <= J; ++j)
    {
        bandlimit = MIN(s2let_bandlimit(j, &parameters), L);
        samples += bandlimit * (2 * bandlimit - 1);
    }

    // Account for the different wavelet normalisations
    for (i = 0; i < samples; ++i)
        f_wav_dir[i] *= sqrt(2*PI);

    // Compute the maximum absolute error in the computed wavelet transform
    wav_error = maxerr_cplx(f_wav_axisym, f_wav_dir, (J-J_min+1)*L*(2*L-1));
    bandlimit = MIN(s2let_bandlimit(J_min-1, &parameters), L);
    scal_error = maxerr_cplx(f_scal_axisym, f_scal_dir, bandlimit*(2*bandlimit-1));

    printf("  - Maximum abs error in wavelets :         %6.5e\n", wav_error);
    printf("  - Maximum abs error in scaling function : %6.5e\n", scal_error);
    fflush(NULL);
}

void s2let_transform_performance_test(int B, int J_min, int NREPEAT, int NSCALE, int seed)
{
  complex double *f, *flm, *flm_rec, *f_rec, *f_wav, *f_scal;
  clock_t time_start, time_end;
  int sc, repeat;
  double tottime_analysis = 0, tottime_synthesis = 0;
  double accuracy = 0.0;

  int L = 2;

  for (sc=0; sc<NSCALE; sc++) {

    L *= 2;

    s2let_lm_allocate(&flm, L);

    printf(" > L =  %i \n", L);
    for (repeat=0; repeat<NREPEAT; repeat++){

      printf("  -> Iteration : %i on %i \n",repeat+1,NREPEAT);

      s2let_lm_random_flm(flm, L, 0, seed);
      s2let_mw_allocate(&f, L);
      s2let_mw_alm2map(f, flm, L);
      s2let_transform_axisym_allocate_mw_f_wav(&f_wav, &f_scal, B, L, J_min);

      time_start = clock();
      s2let_transform_axisym_wav_analysis_mw(f_wav, f_scal, f, B, L, J_min);
      time_end = clock();
      tottime_synthesis += (time_end - time_start) / (double)CLOCKS_PER_SEC;
      //printf("  - Duration for S2LET synthesis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);

      s2let_mw_allocate(&f_rec, L);

      time_start = clock();
      s2let_transform_axisym_wav_synthesis_mw(f_rec, f_wav, f_scal, B, L, J_min);
      time_end = clock();
      tottime_analysis += (time_end - time_start) / (double)CLOCKS_PER_SEC;

      s2let_lm_allocate(&flm_rec, L);
      s2let_mw_map2alm(flm_rec, f_rec, L);

      //printf("  - Duration for S2LET analysis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);
      accuracy += maxerr_cplx(flm, flm_rec, L*L);fflush(NULL);

      free(f);
      free(f_rec);
      free(flm_rec);
      free(f_wav);
      free(f_scal);

    }

    tottime_synthesis = tottime_synthesis / (double)NREPEAT;
    tottime_analysis = tottime_analysis / (double)NREPEAT;
    accuracy = accuracy / (double)NREPEAT;

    printf("  - Average duration for S2LET synthesis  : %5.5f seconds\n", tottime_synthesis);
    printf("  - Average duration for S2LET analysis   : %5.5f seconds\n", tottime_analysis);
    printf("  - Average max error on reconstruction  : %6.5e\n", accuracy);

    free(flm);

  }

}



void s2let_transform_performance_multires_test(int B, int J_min, int NREPEAT, int NSCALE, int seed)
{
  complex double *f, *flm, *flm_rec, *f_rec, *f_wav, *f_scal;
  clock_t time_start, time_end;
  int sc, repeat;
  double tottime_analysis = 0, tottime_synthesis = 0;
  double accuracy = 0.0;

  int L = 2;

  for (sc=0; sc<NSCALE; sc++) {

    L *= 2;

    s2let_lm_allocate(&flm, L);

    printf(" > L =  %i \n", L);
    for (repeat=0; repeat<NREPEAT; repeat++){

      printf("  -> Iteration : %i on %i \n",repeat+1,NREPEAT);

      s2let_lm_random_flm(flm, L, 0, seed);
      s2let_mw_allocate(&f, L);
      s2let_mw_alm2map(f, flm, L);
      s2let_transform_axisym_allocate_mw_f_wav_multires(&f_wav, &f_scal, B, L, J_min);

      time_start = clock();
      s2let_transform_axisym_wav_analysis_mw_multires(f_wav, f_scal, f, B, L, J_min);
      time_end = clock();
      tottime_synthesis += (time_end - time_start) / (double)CLOCKS_PER_SEC;
      //printf("  - Duration for S2LET synthesis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);

      s2let_mw_allocate(&f_rec, L);

      time_start = clock();
      s2let_transform_axisym_wav_synthesis_mw_multires(f_rec, f_wav, f_scal, B, L, J_min);
      time_end = clock();
      tottime_analysis += (time_end - time_start) / (double)CLOCKS_PER_SEC;

      s2let_lm_allocate(&flm_rec, L);
      s2let_mw_map2alm(flm_rec, f_rec, L);

      //printf("  - Duration for S2LET analysis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);
      accuracy += maxerr_cplx(flm, flm_rec, L*L);

      free(f);
      free(f_rec);
      free(flm_rec);
      free(f_wav);
      free(f_scal);

    }

    tottime_synthesis = tottime_synthesis / (double)NREPEAT;
    tottime_analysis = tottime_analysis / (double)NREPEAT;
    accuracy = accuracy / (double)NREPEAT;

    printf("  - Average duration for S2LET multires synthesis  : %5.5f seconds\n", tottime_synthesis);
    printf("  - Average duration for S2LET multires analysis   : %5.5f seconds\n", tottime_analysis);
    printf("  - Average max error on reconstruction  : %6.5e\n", accuracy);

    free(flm);

  }

}




void s2let_transform_lm_performance_test(int B, int J_min, int NREPEAT, int NSCALE, int seed)
{
  complex double *flm, *flm_rec, *f_wav_lm, *f_scal_lm;
  clock_t time_start, time_end;
  int sc, repeat;
  double tottime_analysis = 0, tottime_synthesis = 0;
  double accuracy = 0.0;
  double *wav_lm, *scal_lm;

  int L = 4;

  for (sc=0; sc<NSCALE; sc++) {

    L *= 2;

    s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, B, L);
    s2let_transform_axisym_lm_wav(wav_lm, scal_lm, B, L, J_min);
    s2let_lm_allocate(&flm, L);

    printf(" > L =  %i \n", L);
    for (repeat=0; repeat<NREPEAT; repeat++){

      //printf("  -> Iteration : %i on %i \n",repeat+1,NREPEAT);

      s2let_lm_random_flm(flm, L, 0, seed);

      s2let_transform_axisym_lm_allocate_f_wav(&f_wav_lm, &f_scal_lm, B, L, J_min);

      time_start = clock();
      s2let_transform_axisym_lm_wav_analysis(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, B, L, J_min);
      time_end = clock();
      tottime_synthesis += (time_end - time_start) / (double)CLOCKS_PER_SEC;
      //printf("  - Duration for S2LET synthesis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);

      s2let_lm_allocate(&flm_rec, L);

      time_start = clock();
      s2let_transform_axisym_lm_wav_synthesis(flm_rec, f_wav_lm, f_scal_lm, wav_lm, scal_lm, B, L, J_min);
      time_end = clock();
      tottime_analysis += (time_end - time_start) / (double)CLOCKS_PER_SEC;

      //printf("  - Duration for S2LET analysis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);
      accuracy += maxerr_cplx(flm, flm_rec, L*L);

      free(flm_rec);
      free(f_wav_lm);
      free(f_scal_lm);

    }

    tottime_synthesis = tottime_synthesis / (double)NREPEAT;
    tottime_analysis = tottime_analysis / (double)NREPEAT;
    accuracy = accuracy / (double)NREPEAT;

    //printf("  - Average duration for S2LET synthesis  : %5.2e seconds\n", tottime_synthesis);
    //printf("  - Average duration for S2LET analysis   : %5.2e seconds\n", tottime_analysis);
    printf("  - Average duration for S2LET           : %5.2e seconds\n", (tottime_analysis+tottime_synthesis)/2);
    printf("  - Average max error on reconstruction  : %6.5e\n", accuracy);

    free(flm);
    free(wav_lm);
    free(scal_lm);

  }

}



void s2let_transform_lm_performance_multires_test(int B, int J_min, int NREPEAT, int NSCALE, int seed)
{
  complex double *flm, *flm_rec, *f_wav_lm, *f_scal_lm;
  clock_t time_start, time_end;
  int sc, repeat;
  double tottime_analysis = 0, tottime_synthesis = 0;
  double accuracy = 0.0;
  double *wav_lm, *scal_lm;

  int L = 4;

  for (sc=0; sc<NSCALE; sc++) {

    L *= 2;

    s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, B, L);
    s2let_transform_axisym_lm_wav(wav_lm, scal_lm, B, L, J_min);
    s2let_lm_allocate(&flm, L);

    printf(" > L =  %i \n", L);
    for (repeat=0; repeat<NREPEAT; repeat++){

      //printf("  -> Iteration : %i on %i \n",repeat+1,NREPEAT);

      s2let_lm_random_flm(flm, L, 0, seed);

      s2let_transform_axisym_lm_allocate_f_wav_multires(&f_wav_lm, &f_scal_lm, B, L, J_min);

      time_start = clock();
      s2let_transform_axisym_lm_wav_analysis_multires(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, B, L, J_min);
      time_end = clock();
      tottime_synthesis += (time_end - time_start) / (double)CLOCKS_PER_SEC;
      //printf("  - Duration for S2LET synthesis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);

      s2let_lm_allocate(&flm_rec, L);

      time_start = clock();
      s2let_transform_axisym_lm_wav_synthesis_multires(flm_rec, f_wav_lm, f_scal_lm, wav_lm, scal_lm, B, L, J_min);
      time_end = clock();
      tottime_analysis += (time_end - time_start) / (double)CLOCKS_PER_SEC;

      //printf("  - Duration for S2LET analysis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);
      accuracy += maxerr_cplx(flm, flm_rec, L*L);

      free(flm_rec);
      free(f_wav_lm);
      free(f_scal_lm);

    }

    tottime_synthesis = tottime_synthesis / (double)NREPEAT;
    tottime_analysis = tottime_analysis / (double)NREPEAT;
    accuracy = accuracy / (double)NREPEAT;

    //printf("  - Average duration for S2LET multires synthesis  : %5.2e seconds\n", tottime_synthesis);
    //printf("  - Average duration for S2LET multires analysis   : %5.2e seconds\n", tottime_analysis);
    printf("  - Average duration for S2LET multires  : %5.2e seconds\n", (tottime_analysis+tottime_synthesis)/2);
    printf("  - Average max error on reconstruction  : %6.5e\n", accuracy);

    free(flm);
    free(wav_lm);
    free(scal_lm);

  }

}


int main(int argc, char *argv[])
{
  const int L = 32;
  const int N = 16;
  const int B = 3;
  const int J_min = 0;
  const int spin = 2;

  s2let_parameters_t parameters = {};

  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;
  parameters.N = N;
  parameters.spin = spin;

  // This is too often zero, so we add 1 (zero will result in all random
  // numbers being the same).
  const int seed = (int)((double)clock()/(double)CLOCKS_PER_SEC) + 1;
  int l_min = s2let_el_min(&parameters);

  printf("=====================================================================\n");
  printf("Testing S2LET facilities with the MW sampling\n");
  printf("=====================================================================\n");
  printf("PARAMETERS: ");
  printf("L = %i  N = %i  B = %i  l_wav_min = %i  spin = %i  seed = %i\n",
         L, N, B, l_min, spin, seed);
  //s2let_switch_wavtype(3);
  printf("---------------------------------------------------------------------\n");
  printf("> Testing logfact binomial coefficient implementation...\n");
  // Don't use more than 62 as the argument.
  s2let_binomial_coefficient_test(62);
  printf("=====================================================================\n");
  printf("> Testing axisymmetric kernels...\n");
  s2let_tiling_axisym_test(B, L, J_min);
  printf("---------------------------------------------------------------------\n");
  printf("> Testing directionality components...\n");
  s2let_tiling_direction_test(L, N);
  printf("---------------------------------------------------------------------\n");
  printf("> Testing directional wavelets...\n");
  s2let_tiling_wavelet_test(B, L, J_min, N, spin);
  printf("=====================================================================\n");
  printf("> Testing axisymmetric wavelets in harmonic space...\n");
  s2let_transform_axisym_lm_wav_test(B, L, J_min, seed);
  printf("---------------------------------------------------------------------\n");
  printf("> Testing axisymmetric multiresolution algorithm in harmonic space...\n");
  s2let_transform_axisym_lm_wav_multires_test(B, L, J_min, seed);
  printf("---------------------------------------------------------------------\n");
  printf("> Testing directional wavelets in harmonic space...\n");
  s2let_wav_transform_harmonic_test(B, L, J_min, N, spin, seed);
  printf("---------------------------------------------------------------------\n");
  printf("> Testing directional multiresolution algorithm in harmonic space...\n");
  s2let_wav_transform_harmonic_multires_test(B, L, J_min, N, spin, seed);
  printf("=====================================================================\n");
  printf("> Testing axisymmetric wavelets in pixel space...\n");
  s2let_transform_axisym_wav_test(B, L, J_min, seed);
  printf("---------------------------------------------------------------------\n");
  printf("> Testing axisymmetric multiresolution algorithm in pixel space...\n");
  s2let_transform_axisym_wav_multires_test(B, L, J_min, seed);
  printf("---------------------------------------------------------------------\n");
  printf("> Testing directional wavelets in pixel space...\n");
  s2let_wav_transform_mw_test(B, L, J_min, N, spin, seed);
  printf("---------------------------------------------------------------------\n");
  printf("> Testing directional multiresolution algorithm in pixel space...\n");
  s2let_wav_transform_mw_multires_test(B, L, J_min, N, spin, seed);
  printf("=====================================================================\n");
  printf("> Comparing directional and axisymmetric algorithm in pixel space...\n");
  s2let_transform_axisym_vs_directional_mw_test(B, L, J_min, seed);
  printf("---------------------------------------------------------------------\n");
  printf("> Comparing directional and axisymmetric multiresolution algorithm\n");
  printf("  in pixel space...\n");
  s2let_transform_axisym_vs_directional_mw_multires_test(B, L, J_min, seed);
  printf("=====================================================================\n");
  printf("> Testing real axisymmetric wavelets in pixel space...\n");
  s2let_transform_axisym_wav_real_test(B, L, J_min, seed);
  printf("---------------------------------------------------------------------\n");
  printf("> Testing axisymmetric multiresolution algorithm for real function...\n");
  s2let_transform_axisym_wav_multires_real_test(B, L, J_min, seed);
  printf("---------------------------------------------------------------------\n");
  printf("> Testing directional wavelets in pixel space for real function...\n");
  s2let_wav_transform_mw_real_test(B, L, J_min, N, seed);
  printf("---------------------------------------------------------------------\n");
  printf("> Testing directional multiresolution algorithm for real function...\n");
  s2let_wav_transform_mw_multires_real_test(B, L, J_min, N, seed);

    /*
  const int NREPEAT = 50;
  const int NSCALE = 9;
    printf("==========================================================\n");
    printf("> Testing performances at full resolution...\n");
    s2let_transform_lm_performance_test(B, J_min, NREPEAT, NSCALE, seed);
    printf("----------------------------------------------------------\n");
    printf("> Testing performances with multiresolution...\n");
    s2let_transform_lm_performance_multires_test(B, J_min, NREPEAT, NSCALE, seed);
    printf("----------------------------------------------------------\n");
    printf("> Testing performances at full resolution...\n");
    //s2let_transform_performance_test(B, J_min, NREPEAT, NSCALE, seed);
    printf("----------------------------------------------------------\n");
    printf("> Testing performances with multiresolution...\n");
    //s2let_transform_performance_multires_test(B, J_min, NREPEAT, NSCALE, seed);
    */
  printf("==========================================================\n");



  return 0;
}

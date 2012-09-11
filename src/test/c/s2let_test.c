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
void s2let_axisym_tiling_test(int B, int L, int J_min)
{		
	double *kappa, *kappa0;

	// Allocate the kernels corresponding to the parameters B, L
	s2let_axisym_allocate_tiling(&kappa, &kappa0, B, L);

	// Construct the tiling of harmonic space
	s2let_axisym_tiling(kappa, kappa0, B, L, J_min);

	// Check that they recover the identity relation,
	// ensuring exactness of the wavelet transform.
	double res = s2let_axisym_check_identity(kappa, kappa0, B, L, J_min);
	printf("  - Identity residuals : %6.5e\n", res);

	free(kappa);
	free(kappa0);
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
void s2let_axisym_wav_lm_test(int B, int L, int J_min, int seed)
{
	clock_t time_start, time_end;
	double *wav_lm, *scal_lm;

	// Allocate the wavelet kernels
	s2let_axisym_allocate_wav_lm(&wav_lm, &scal_lm, B, L);

	// Compute the wavelet kernels
	time_start = clock();
	s2let_axisym_wav_lm(wav_lm, scal_lm, B, L, J_min);
	time_end = clock();
	printf("  - Generate wavelets  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);
	 
	complex double *f_wav_lm, *f_scal_lm, *flm, *flm_rec;
	s2let_allocate_lm(&flm, L);
	s2let_allocate_lm(&flm_rec, L);

	// Generate a random spherical harmonic decomposition
	s2let_axisym_random_flm(flm, L, seed);
	
	// Allocate space for the wavelet scales (their harmonic coefficients)
	s2let_axisym_allocate_f_wav_lm(&f_wav_lm, &f_scal_lm, B, L, J_min);
	
	// Perform the wavelet transform through exact harmonic tiling
	time_start = clock();
	s2let_axisym_wav_analysis_lm(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	// Reconstruct the initial harmonic coefficients from those of the wavelets
	time_start = clock();
	s2let_axisym_wav_synthesis_lm(flm_rec, f_wav_lm, f_scal_lm, wav_lm, scal_lm, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	// Compute the maximum absolute error on the harmonic coefficients
	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flm, flm_rec, L*L));

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
void s2let_axisym_wav_lm_multires_test(int B, int L, int J_min, int seed)
{
	clock_t time_start, time_end;
	double *wav_lm, *scal_lm;

	// Allocate the wavelet kernels
	s2let_axisym_allocate_wav_lm(&wav_lm, &scal_lm, B, L);

	// Compute the wavelet kernels
	time_start = clock();
	s2let_axisym_wav_lm(wav_lm, scal_lm, B, L, J_min);
	time_end = clock();
	printf("  - Generate wavelets  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);
	 
	complex double *f_wav_lm, *f_scal_lm, *flm, *flm_rec;
	s2let_allocate_lm(&flm, L);
	s2let_allocate_lm(&flm_rec, L);

	// Generate a random spherical harmonic decomposition
	s2let_axisym_random_flm(flm, L, seed);
	
	// Allocate space for the wavelet scales (their harmonic coefficients)
	s2let_axisym_allocate_f_wav_multires_lm(&f_wav_lm, &f_scal_lm, B, L, J_min);
	
	// Perform the wavelet transform through exact harmonic tiling
	time_start = clock();
	s2let_axisym_wav_analysis_multires_lm(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	// Reconstruct the initial harmonic coefficients from those of the wavelets
	time_start = clock();
	s2let_axisym_wav_synthesis_multires_lm(flm_rec, f_wav_lm, f_scal_lm, wav_lm, scal_lm, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);
	
	// Compute the maximum absolute error on the harmonic coefficients
	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flm, flm_rec, L*L));

	free(flm);
	free(flm_rec);
	free(f_wav_lm);
	free(f_scal_lm);
	free(wav_lm);
	free(scal_lm);
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
void s2let_axisym_wav_test(int B, int L, int J_min, int seed)
{
	clock_t time_start, time_end;
	int spin = 0;
	int verbosity = 0;
	ssht_dl_method_t dl_method = SSHT_DL_RISBO;
	//int J = s2let_j_max(L, B);

	complex double *f, *f_rec, *flm, *flm_rec;
	s2let_allocate_lm(&flm, L);
	s2let_allocate_lm(&flm_rec, L);
	s2let_allocate_mw(&f, L);
	s2let_allocate_mw(&f_rec, L);

	// Generate random harmonic coefficients for a complex signal
	s2let_axisym_random_flm(flm, L, seed);

	// Construct the corresponding signal on the sphere (MW sampling)
 	ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

 	// Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
	complex double *f_wav, *f_scal;
	s2let_axisym_allocate_f_wav(&f_wav, &f_scal, B, L, J_min);

	// Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
	time_start = clock();
	s2let_axisym_wav_analysis(f_wav, f_scal, f, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	// Reconstruct the initial signal from the wavelet maps from scratch
	time_start = clock();
	s2let_axisym_wav_synthesis(f_rec, f_wav, f_scal, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	// Compute the initial harmonic coefficients back
	ssht_core_mw_forward_sov_conv_sym(flm_rec, f_rec, L, spin, dl_method, verbosity);

	// Compute the maximum absolute error on the harmonic coefficients
	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flm, flm_rec, L*L));

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
void s2let_axisym_wav_real_test(int B, int L, int J_min, int seed)
{
	clock_t time_start, time_end;
	int verbosity = 0;
	ssht_dl_method_t dl_method = SSHT_DL_RISBO;
	//int J = s2let_j_max(L, B);

	complex *flm, *flm_rec;
	double *f, *f_rec;
	s2let_allocate_lm(&flm, L);
	s2let_allocate_lm(&flm_rec, L);
	s2let_allocate_mw_real(&f, L);
	s2let_allocate_mw_real(&f_rec, L);

	// Generate random harmonic coefficients for a real signal
	s2let_axisym_random_flm_real(flm, L, seed);

	// Construct the corresponding signal on the sphere (MW sampling)
 	ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

	// Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
	double *f_wav, *f_scal;
	s2let_axisym_allocate_f_wav_real(&f_wav, &f_scal, B, L, J_min);

	// Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
	time_start = clock();
	s2let_axisym_wav_analysis_real(f_wav, f_scal, f, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	// Reconstruct the initial signal from the wavelet maps from scratch
	time_start = clock();
	s2let_axisym_wav_synthesis_real(f_rec, f_wav, f_scal, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	// Compute the initial harmonic coefficients back
	ssht_core_mw_forward_sov_conv_sym_real(flm_rec, f_rec, L, dl_method, verbosity);

	// Compute the maximum absolute error on the harmonic coefficients
	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flm, flm_rec, L*L));
	
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
void s2let_axisym_wav_multires_test(int B, int L, int J_min, int seed)
{
	clock_t time_start, time_end;
	int spin = 0;
	int verbosity = 0;
	ssht_dl_method_t dl_method = SSHT_DL_RISBO;
	//int J = s2let_j_max(L, B);

	complex double *f, *f_rec, *flm, *flm_rec;
	s2let_allocate_lm(&flm, L);
	s2let_allocate_lm(&flm_rec, L);
	s2let_allocate_mw(&f, L);
	s2let_allocate_mw(&f_rec, L);

	// Generate random harmonic coefficients for a complex signal
	s2let_axisym_random_flm(flm, L, seed);

	// Construct the corresponding signal on the sphere (MW sampling)
 	ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

	// Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
	complex double *f_wav, *f_scal;
	s2let_axisym_allocate_f_wav_multires(&f_wav, &f_scal, B, L, J_min);
	
	// Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
	time_start = clock();
	s2let_axisym_wav_analysis_multires(f_wav, f_scal, f, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	// Reconstruct the initial signal from the wavelet maps from scratch
	time_start = clock();
	s2let_axisym_wav_synthesis_multires(f_rec, f_wav, f_scal, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	// Compute the initial harmonic coefficients back
	ssht_core_mw_forward_sov_conv_sym(flm_rec, f_rec, L, spin, dl_method, verbosity);

	// Compute the maximum absolute error on the harmonic coefficients
	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flm, flm_rec, L*L));
	
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
void s2let_axisym_wav_multires_real_test(int B, int L, int J_min, int seed)
{
	clock_t time_start, time_end;
	int verbosity = 0;
	ssht_dl_method_t dl_method = SSHT_DL_RISBO;
	//int J = s2let_j_max(L, B);

	complex *flm, *flm_rec;
	double *f, *f_rec;
	s2let_allocate_lm(&flm, L);
	s2let_allocate_lm(&flm_rec, L);
	s2let_allocate_mw_real(&f, L);
	s2let_allocate_mw_real(&f_rec, L);

	// Generate random harmonic coefficients for a real signal
	s2let_axisym_random_flm_real(flm, L, seed);

	// Construct the corresponding signal on the sphere (MW sampling)
 	ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

	// Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
	double *f_wav, *f_scal;
	s2let_axisym_allocate_f_wav_multires_real(&f_wav, &f_scal, B, L, J_min);

	// Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
	time_start = clock();
	s2let_axisym_wav_analysis_multires_real(f_wav, f_scal, f, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	// Reconstruct the initial signal from the wavelet maps from scratch
	time_start = clock();
	s2let_axisym_wav_synthesis_multires_real(f_rec, f_wav, f_scal, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	// Compute the initial harmonic coefficients back
	ssht_core_mw_forward_sov_conv_sym_real(flm_rec, f_rec, L, dl_method, verbosity);

	// Compute the maximum absolute error on the harmonic coefficients
	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flm, flm_rec, L*L));

	free(f);
	free(f_rec);
	free(f_wav);
	free(f_scal);
}


int main(int argc, char *argv[]) 
{
	
	const int L = 64;
	const int B = 3;
	const int J_min = 1;
	const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);
	int l_min = s2let_el_min(B, J_min);

	printf("==========================================================\n");
	printf("Testing S2LET facilities with the MW sampling\n");
	printf("==========================================================\n");
	printf("PARAMETERS: ");
	printf("  L = %i   B = %i   l_wav_min = %i   seed = %i\n", L, B, l_min, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing harmonic tiling...\n");
	s2let_axisym_tiling_test(B, L, J_min);
	printf("==========================================================\n");
	printf("> Testing axisymmetric wavelets in harmonics space...\n");
	s2let_axisym_wav_lm_test(B, L, J_min, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing multiresolution algorithm in harmonics space...\n");
	s2let_axisym_wav_lm_multires_test(B, L, J_min, seed);
	printf("==========================================================\n");
	printf("> Testing axisymmetric wavelets in pixel space...\n");
	s2let_axisym_wav_test(B, L, J_min, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing multiresolution algorithm...\n");
	s2let_axisym_wav_multires_test(B, L, J_min, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing real axisymmetric wavelets in pixel space...\n");
	s2let_axisym_wav_real_test(B, L, J_min, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing multiresolution algorithm for real function...\n");
	s2let_axisym_wav_multires_real_test(B, L, J_min, seed);
	printf("==========================================================\n");


	return 0;		
}

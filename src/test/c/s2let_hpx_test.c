// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <assert.h>
#include <complex.h> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <fftw3.h> 
#include <ssht.h>

/*!
 * Perform HEALPIX spherical harmonic transform back and forth.
 * Test that the interfaces to HEALPIX routines work.
 *
 * \param[out]  nside Healpix resolution.
 * \param[in]  L Band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_axisym_hpx_test(int nside, int L, int seed)
{
	double *f, *f_rec;
	complex double *flm, *flm_rec;
	s2let_allocate_lm(&flm, L);
	s2let_allocate_lm(&flm_rec, L);
	s2let_allocate_hpx_real(&f, nside);
	s2let_allocate_hpx_real(&f_rec, nside);

	// Generate random harmonic coefficients
	s2let_axisym_random_flm_real(flm, L, seed);

	// Reconstruct the corresponding signal on the sphere on a healpix map
	s2let_hpx_alm2map_real(f, flm, nside, L);

	// Decompose it again to get the harmonic coefficients back
	s2let_hpx_map2alm_real(flm_rec, f, nside, L);

	// Compute the maximum absolute error on the harmonic coefficients
	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flm, flm_rec, L*L));
	
	free(f);
	free(f_rec);
	free(flm);
	free(flm_rec);
}


/*!
 * Perform wavelet transform with HEALPIX map back and forth.
 * Test that the interfaces to HEALPIX routines work.
 *
 * \param[out]  nside Healpix resolution.
 * \param[in]  L Band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_axisym_hpx_wav_test(int nside, int B, int L, int J_min, int seed)
{
	clock_t time_start, time_end;

	double *f, *f_rec;
	complex double *flm, *flm_rec;
	s2let_allocate_lm(&flm, L);
	s2let_allocate_lm(&flm_rec, L);
	s2let_allocate_hpx_real(&f, nside);
	s2let_allocate_hpx_real(&f_rec, nside);

	// Generate random harmonic coefficients
	s2let_axisym_random_flm_real(flm, L, seed);

	// Reconstruct the corresponding signal on the sphere on a healpix map
	s2let_hpx_alm2map_real(f, flm, nside, L);

	// Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
	double *f_wav, *f_scal;
	s2let_axisym_hpx_allocate_f_wav_real(&f_wav, &f_scal, nside, B, L, J_min);

	// Perform wavelet analysis from scratch with all signals given on the sphere (Healpix sampling)
	time_start = clock();
	s2let_axisym_hpx_wav_analysis_real(f_wav, f_scal, f, nside, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	// Reconstruct the initial healpix map from the wavelet healpix maps
	time_start = clock();
	s2let_axisym_hpx_wav_synthesis_real(f_rec, f_wav, f_scal, nside, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	// Get the harmonic coefficients back
	s2let_hpx_map2alm_real(flm_rec, f_rec, nside, L);

	// Compute the maximum absolute error on the harmonic coefficients
	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flm, flm_rec, L*L));
	
	free(f);
	free(f_rec);
	free(f_wav);
	free(f_scal);
}

/*!
 * Test the Input-Output facilities for Healpix maps (to FITS filts)
 *
 * \param[in]  nside Healpix resolution.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_hpx_io_test(int nside, int L, int seed)
{
	double *f, *f_rec;
	complex double *flm, *flm_rec;
	s2let_allocate_lm(&flm, L);
	s2let_allocate_lm(&flm_rec, L);
	s2let_allocate_hpx_real(&f, nside);
	s2let_allocate_hpx_real(&f_rec, nside);

	// Generate random harmonic coefficients
	s2let_axisym_random_flm_real(flm, L, seed);

	// Construct the corresponding real signal on a healpix map
	s2let_hpx_alm2map_real(f, flm, nside, L);

	char file[100] = "temp.fits";

	// Remove the file if it exists
	remove(file);
	// Write the signal to file
	s2let_write_hpx_map(file, f, nside);

	// Read the signal from file
	s2let_read_hpx_map(f_rec, file, nside);
	// Clean
	remove(file);

	// Get the harmonic coefficients back 
	s2let_hpx_map2alm_real(flm_rec, f, nside, L);

	// Compute the maximum absolute error on the harmonic coefficients
	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flm, flm_rec, L*L));
	
	free(f);
	free(f_rec);
	free(flm);
	free(flm_rec);
}

/*!
 * Test the Input-Output facilities for the MW sampling (to FITS filts)
 *
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_mw_io_test(int L, int seed)
{
	int verbosity = 0;
	ssht_dl_method_t dl_method = SSHT_DL_RISBO;

	double *f, *f_rec;
	complex double *flm, *flm_rec;
	s2let_allocate_lm(&flm, L);
	s2let_allocate_lm(&flm_rec, L);
	s2let_allocate_mw_real(&f, L);
	s2let_allocate_mw_real(&f_rec, L);

	// Generate random harmonic coefficients
	s2let_axisym_random_flm_real(flm, L, seed);

	// Construct the corresponding real signal, on MW sampling 
	ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

	char file[100] = "temp.fits";

	// Remove the file if it exists
	remove(file);
	// Write the signal to file
	s2let_write_mw_map(file, f, L);

	// Read the band-limit from file
	int Lread = s2let_read_mw_bandlimit(file);

	// Read the signal from file
	s2let_read_mw_map(f_rec, file, Lread);
	// Clean
	remove(file);

	// Get the harmonic coefficients back 
	ssht_core_mw_forward_sov_conv_sym_real(flm_rec, f_rec, L, dl_method, verbosity);

	// Compute the maximum absolute error on the harmonic coefficients
	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flm, flm_rec, L*L));
	
	free(f);
	free(f_rec);
	free(flm);
	free(flm_rec);
}

int main(int argc, char *argv[]) 
{
	
	const int L = 256;
	const int B = 2;
	const int J_min = 2;
	int nside = 128;
	const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);

	printf("==========================================================\n");
	printf("Testing S2LET facilities with the HEALPIX sampling \n");
	printf("==========================================================\n");
	printf("PARAMETERS: ");
	printf("  L = %i   B = %i   nside = %i   seed = %i\n", L, B, nside, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing interface to healpix transform...\n");
	s2let_axisym_hpx_test(nside, L, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing real axisymmetric wavelets in pixel space...\n");
	s2let_axisym_hpx_wav_test(nside, B, L, J_min, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing IO functions for MW sampling...\n");
	s2let_mw_io_test(L, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing IO functions for HEALPIX sampling...\n");
	s2let_hpx_io_test(nside, L, seed);
	printf("==========================================================\n");


	return 0;		
}

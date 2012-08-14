// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <assert.h>

/*!
 * Max absolute error between two complex arrays
 */
double maxerr_cplx(complex double *a, complex double *b, int size)
{
	double value = 0;
	int i;
	for(i = 0; i<size; i++){
		value = MAX( cabs( a[i]-b[i] ), value );
	}
	return value;
}

/*!
 * Max absolute error between two real arrays
 */
double maxerr(double *a, double *b, int size)
{
	double value = 0;
	int i;
	for(i = 0; i<size; i++){
		value = MAX( abs( a[i]-b[i] ), value );
	}
	return value;
}

/*!
 * Random number from seed (Numerical Recipes)
 */
double ran2_dp(int idum) {
  int IM1=2147483563,IM2=2147483399,IMM1=IM1-1, 
    IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, 
    NTAB=32,NDIV=1+IMM1/NTAB;

  double AM=1./IM1,EPS=1.2e-7,RNMX=1.-EPS;
  int j,k;
  static int iv[32],iy,idum2 = 123456789; 
  // N.B. in C static variables are initialised to 0 by default.

  if (idum <= 0) {
    idum= (-idum>1 ? -idum : 1); // max(-idum,1);
    idum2=idum;
    for(j=NTAB+8;j>=1;j--) {
      k=idum/IQ1;
      idum=IA1*(idum-k*IQ1)-k*IR1;
      if (idum < 0) idum=idum+IM1;
      if (j < NTAB) iv[j-1]=idum;
    }
    iy=iv[0];
  }
  k=idum/IQ1;
  idum=IA1*(idum-k*IQ1)-k*IR1;
  if (idum < 0) idum=idum+IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2=idum2+IM2;
  j=1+iy/NDIV;
  iy=iv[j-1]-idum2;
  iv[j-1]=idum;
  if(iy < 1)iy=iy+IMM1;
  return (AM*iy < RNMX ? AM*iy : RNMX); // min(AM*iy,RNMX);
}

/*!
 * Generate random harmonic coefficients for a complex map.
 *
 * \param[out]  flm Harmonic coefficients.
 * \param[in]  L Band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_axisym_random_flm(complex double *flm, int L, int seed)
{
	int i;
	srand( time(NULL) );
	for (i=0; i<L*L; i++){
		flm[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
	}
}

/*!
 * Generate random harmonic coefficients corresponding to a real map.
 *
 * \param[out]  flm Harmonic coefficients.
 * \param[in]  L Band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_axisym_random_flm_real(complex double *flm, int L, int seed) {

  int el, m, msign, i, i_op;
  for (el=0; el<L; el++) {
    m = 0;
    i = el*el + el + m ;
    flm[i] = (2.0*ran2_dp(seed) - 1.0);
    for (m=1; m<=el; m++) {
      i = el*el + el + m ;
      flm[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
      i_op = el*el + el - m ;
      msign = m & 1;
      msign = 1 - msign - msign; // (-1)^m
      flm[i_op] = msign * conj(flm[i]);
    }
  }
}

/*!
 * Test the identity relation of the wavelet tilling in harmonic space.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_tilling_test(int B, int L, int J_min)
{		
	double *kappa, *kappa0;

	// Allocate the kernels corresponding to the parameters B, L
	s2let_axisym_allocate_tilling(&kappa, &kappa0, B, L);

	// Construct the tilling of harmonic space
	s2let_axisym_tilling(kappa, kappa0, B, L, J_min);

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
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	flm_rec = (complex double*)calloc(L * L, sizeof(complex double));

	// Generate a random spherical harmonic decomposition
	s2let_axisym_random_flm(flm, L, seed);
	
	// Allocate space for the wavelet scales (their harmonic coefficients)
	s2let_axisym_allocate_f_wav_lm(&f_wav_lm, &f_scal_lm, B, L, J_min);
	
	// Perform the wavelet transform through exact harmonic tilling
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
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	flm_rec = (complex double*)calloc(L * L, sizeof(complex double));

	// Generate a random spherical harmonic decomposition
	s2let_axisym_random_flm(flm, L, seed);
	
	// Allocate space for the wavelet scales (their harmonic coefficients)
	s2let_axisym_allocate_f_wav_multires_lm(&f_wav_lm, &f_scal_lm, B, L, J_min);
	
	// Perform the wavelet transform through exact harmonic tilling
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
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	flm_rec = (complex double*)calloc(L * L, sizeof(complex double));
	f = (complex double*)calloc(L*(2*L-1), sizeof(complex double));
	f_rec = (complex double*)calloc(L*(2*L-1), sizeof(complex double));

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
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	flm_rec = (complex double*)calloc(L * L, sizeof(complex double));
	f = (double*)calloc(L*(2*L-1), sizeof(double));
	f_rec = (double*)calloc(L*(2*L-1), sizeof(double));

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
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	flm_rec = (complex double*)calloc(L * L, sizeof(complex double));
	f = (complex double*)calloc(L*(2*L-1), sizeof(complex double));
	f_rec = (complex double*)calloc(L*(2*L-1), sizeof(complex double));

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
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	flm_rec = (complex double*)calloc(L * L, sizeof(complex double));
	f = (double*)calloc(L*(2*L-1), sizeof(double));
	f_rec = (double*)calloc(L*(2*L-1), sizeof(double));

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
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	f = (double*)calloc(L * (2*L-1), sizeof(double));
	f_rec = (double*)calloc(L * (2*L-1), sizeof(double));

	// Generate random harmonic coefficients
	s2let_axisym_random_flm_real(flm, L, seed);

	// Construct the corresponding real signal, on MW sampling 
	ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

	char file[100] = "temp.fits";

	// Remove the file if it exists
	remove(file);
	// Write the signal to file
	write_mw_map(file, f, L);

	// Read the band-limit from file
	int Lread = read_mw_bandlimit(file);

	// Read the signal from file
	read_mw_map(f_rec, file, Lread);
	// Clean
	remove(file);

	// Get the harmonic coefficients back 
	flm_rec = (complex double*)calloc(L * L, sizeof(complex double));
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
	
	const int L = 64;
	const int B = 2;
	const int J_min = 0;
	const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);
	int l_min = s2let_el_min(B, J_min);

	printf("==========================================================\n");
	printf("Testing S2LET facilities with the MW sampling\n");
	printf("==========================================================\n");
	printf("PARAMETERS: ");
	printf("  L = %i   B = %i   l_wav_min = %i   seed = %i\n", L, B, l_min, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing IO functions for MW sampling...\n");
	s2let_mw_io_test(L, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing harmonic tilling...\n");
	s2let_axisym_tilling_test(B, L, J_min);
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

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
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	flm_rec = (complex double*)calloc(L * L, sizeof(complex double));
	f = (double*)calloc(12*nside*nside, sizeof(double));
	f_rec = (double*)calloc(12*nside*nside, sizeof(double));

	// Generate random harmonic coefficients
	s2let_axisym_random_flm_real(flm, L, seed);

	// Reconstruct the corresponding signal on the sphere on a healpix map
	healpix_inverse_real(f, flm, nside, L);

	// Decompose it again to get the harmonic coefficients back
	healpix_forward_real(flm_rec, f, nside, L);

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
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	flm_rec = (complex double*)calloc(L * L, sizeof(complex double));
	f = (double*)calloc(12*nside*nside, sizeof(double));
	f_rec = (double*)calloc(12*nside*nside, sizeof(double));

	// Generate random harmonic coefficients
	s2let_axisym_random_flm_real(flm, L, seed);

	// Reconstruct the corresponding signal on the sphere on a healpix map
	healpix_inverse_real(f, flm, nside, L);

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
	healpix_forward_real(flm_rec, f_rec, nside, L);

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
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	flm_rec = (complex double*)calloc(L * L, sizeof(complex double));
	f = (double*)calloc(12*nside*nside, sizeof(double));
	f_rec = (double*)calloc(12*nside*nside, sizeof(double));

	// Generate random harmonic coefficients
	s2let_axisym_random_flm_real(flm, L, seed);

	// Construct the corresponding real signal on a healpix map
	healpix_inverse_real(f, flm, nside, L);

	char file[100] = "temp.fits";

	// Remove the file if it exists
	remove(file);
	// Write the signal to file
	write_healpix_map(file, f, nside);

	// Read the signal from file
	read_healpix_map(f_rec, file, nside);
	// Clean
	remove(file);

	// Get the harmonic coefficients back 
	healpix_forward_real(flm_rec, f, nside, L);

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
	const int J_min = 2;
	int nside = 32;
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
	printf("> Testing interface to IO-fits healpix functions...\n");
	s2let_hpx_io_test(nside, L, seed);
	printf("==========================================================\n");


	return 0;		
}

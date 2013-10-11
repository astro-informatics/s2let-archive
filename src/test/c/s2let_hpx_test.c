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
void s2let_axisym_hpx_test(double *accuracy, double *timing, int nside, int L, int seed)
{
  clock_t time_start, time_end;
  double *f, *f_rec;
  complex double *flm, *flm_rec;
  s2let_lm_allocate(&flm, L);
  s2let_lm_allocate(&flm_rec, L);
  s2let_hpx_allocate_real(&f, nside);
  s2let_hpx_allocate_real(&f_rec, nside);

  // Generate random harmonic coefficients
  s2let_lm_random_flm_real(flm, L, seed);

  fflush(NULL);time_start = clock();
  // Reconstruct the corresponding signal on the sphere on a healpix map
  s2let_hpx_alm2map_real(f, flm, nside, L);
  // Decompose it again to get the harmonic coefficients back
  s2let_hpx_map2alm_real(flm_rec, f, nside, L);
  time_end = clock();fflush(NULL);

  int l, m;
  for(l = 0; l<L; l++){
  for(m = 0; m<=l; m++){
      //printf("l=%i, m=%i, val=%f \n",l,m,cabs(flm_rec[l*l+l+m]-flm[l*l+l+m]));
    }}


  *timing += (time_end - time_start) / 2 / (double)CLOCKS_PER_SEC ;
  *accuracy = maxerr_cplx(flm, flm_rec, L*L);

  // Compute the maximum absolute error on the harmonic coefficients
  //printf("  - Maximum abs error  : %6.5e\n", accuracy);
	
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
void s2let_axisym_hpx_wav_test(double *accuracy, double *timing, int nside, int B, int L, int J_min, int seed)
{
  clock_t time_start, time_end;

  double *f, *f_rec;
  complex double *flm, *flm_rec;
  s2let_lm_allocate(&flm, L);
  s2let_lm_allocate(&flm_rec, L);
  s2let_hpx_allocate_real(&f, nside);
  s2let_hpx_allocate_real(&f_rec, nside);

  // Generate random harmonic coefficients
  s2let_lm_random_flm_real(flm, L, seed);

  // Reconstruct the corresponding signal on the sphere on a healpix map
  s2let_hpx_alm2map_real(f, flm, nside, L);

  // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
  double *f_wav, *f_scal;
  s2let_axisym_hpx_allocate_f_wav_real(&f_wav, &f_scal, nside, B, L, J_min);

  // Perform wavelet analysis from scratch with all signals given on the sphere (Healpix sampling)
  fflush(NULL);time_start = clock();
  s2let_axisym_hpx_wav_analysis_real(f_wav, f_scal, f, nside, B, L, J_min);
  time_end = clock();fflush(NULL);
  //printf("  - Wavelet analysis   : %4.4f seconds\n", 
	// (time_end - time_start) / (double)CLOCKS_PER_SEC);
  *timing = (time_end - time_start) / 2 / (double)CLOCKS_PER_SEC ;

  // Reconstruct the initial healpix map from the wavelet healpix maps
  time_start = clock();
  s2let_axisym_hpx_wav_synthesis_real(f_rec, f_wav, f_scal, nside, B, L, J_min);
  time_end = clock();
  //printf("  - Wavelet synthesis  : %4.4f seconds\n", 
	// (time_end - time_start) / (double)CLOCKS_PER_SEC);
  *timing += (time_end - time_start) / 2 / (double)CLOCKS_PER_SEC ;

  // Get the harmonic coefficients back
  s2let_hpx_map2alm_real(flm_rec, f_rec, nside, L);

  // Compute the maximum absolute error on the harmonic coefficients
  // printf("  - Maximum abs error  : %6.5e\n", 
	// maxerr_cplx(flm, flm_rec, L*L));

  *accuracy = maxerr_cplx(flm, flm_rec, L*L);
	
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
  s2let_lm_allocate(&flm, L);
  s2let_lm_allocate(&flm_rec, L);
  s2let_hpx_allocate_real(&f, nside);
  s2let_hpx_allocate_real(&f_rec, nside);

  // Generate random harmonic coefficients
  s2let_lm_random_flm_real(flm, L, seed);

  // Construct the corresponding real signal on a healpix map
  s2let_hpx_alm2map_real(f, flm, nside, L);

  char file[100] = "temp.fits";

  // Remove the file if it exists
  remove(file);
  // Write the signal to file
  s2let_hpx_write_map(file, f, nside);

  // Read the signal from file
  s2let_hpx_read_map(f_rec, file, nside);
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
  s2let_lm_allocate(&flm, L);
  s2let_lm_allocate(&flm_rec, L);
  s2let_mw_allocate_real(&f, L);
  s2let_mw_allocate_real(&f_rec, L);

  // Generate random harmonic coefficients
  s2let_lm_random_flm_real(flm, L, seed);

  // Construct the corresponding real signal, on MW sampling 
  ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

  char file[100] = "temp.fits";

  // Remove the file if it exists
  remove(file);
  // Write the signal to file
  s2let_fits_mw_write_map(file, f, L);

  // Read the band-limit from file
  int Lread = s2let_fits_mw_read_bandlimit(file);

  // Read the signal from file
  s2let_fits_mw_read_map(f_rec, file, Lread);
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
	
  int L = 8;
  const int B = 2;
  const int J_min = 0;
  int J = s2let_j_max(L, B);
  int nside = 32;
  int repeat, NREPEAT = 10;
  double timing, accuracy, timing_tot, accuracy_tot;
  const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);

  printf("==========================================================\n");
  printf("Testing S2LET facilities with the HEALPIX sampling \n");
  printf("==========================================================\n");
  printf("PARAMETERS: ");
  printf("  L = %i   B = %i   nside = %i   seed = %i\n", L, B, nside, seed);
  printf("----------------------------------------------------------\n");
  printf("> Testing IO functions for MW sampling...\n");
  s2let_mw_io_test(L, seed);
  printf("----------------------------------------------------------\n");
  printf("> Testing IO functions for HEALPIX sampling...\n");
  s2let_hpx_io_test(nside, L, seed);
  printf("==========================================================\n");


  printf("> Extensive test of accuracy for various Nside and lmax..\n");
  const int i_nmin = 4;
  const int i_nmax = 7;
  int n, l, i_lmin, i_lmax;

  /*FILE *file1, *file2;
  file1 = fopen("s2let_hpx_full.txt", "w");
  if (file1 == NULL) {
         printf("I couldn't open s2let_hpx.txt for writing.\n");
         exit(0);
      }
  */

  for (n=i_nmin; n<=i_nmax; n++){
  for (l=-1; l<=2; l++){  
    nside = pow(2, n);
    L = pow(2,n+l);
    if (l == 2) L = 3*nside;
      accuracy_tot = 0.0;
      timing_tot = 0.0;
      for (repeat=0; repeat<NREPEAT; repeat++){
        s2let_axisym_hpx_test(&accuracy, &timing, nside, L, seed);
        accuracy_tot += accuracy;
        timing_tot += timing;
      }
      accuracy_tot /= NREPEAT;
      timing_tot /= NREPEAT;
      printf("1-HPX: Nside = %i ; Lmax = %i : Accuracy : %6.3e in %6.3e sec \n",nside, L, accuracy_tot,timing_tot);
      //fprintf(file1, " %i %i %6.3e %6.3e \n", nside, L, accuracy_tot,timing_tot);

    }
  } 
  //close(file1);

  printf("==========================================================\n");

/*
  const int i_nmin = 5;
  const int i_nmax = 10;
  int n, l, i_lmin, i_lmax;

  FILE *file1, *file2;
  file1 = fopen("s2let_hpx.txt", "w");
  if (file1 == NULL) {
         printf("I couldn't open s2let_hpx.txt for writing.\n");
         exit(0);
      }
  file2 = fopen("s2let_hpx_wav.txt", "w");
  if (file2 == NULL) {
         printf("I couldn't open s2let_hpx_wav.txt for writing.\n");
         exit(0);
      }

  printf("> Testing real axisymmetric wavelets in pixel space...\n");
  
   printf("NbrScale = %i\n",J);
  for (n=i_nmin; n<=i_nmax; n++){

    if (n > 9) NREPEAT = 4;
    nside = pow(2, n);
    i_lmin = n-2;
    i_lmax = n+1; 
    for (l=i_lmin; l<=i_lmax; l++){
      L = pow(2,l);

      accuracy_tot = 0.0;
      timing_tot = 0.0;
      for (repeat=0; repeat<NREPEAT; repeat++){
        s2let_axisym_hpx_test(&accuracy, &timing, nside, L, seed);
        accuracy_tot += accuracy;
        timing_tot += timing;
      }
      accuracy_tot /= NREPEAT;
      timing_tot /= NREPEAT;
      printf("1-HPX: Nside = %i ; Lmax = %i : Accuracy : %6.3e in %6.3e sec \n",nside, L, accuracy_tot,timing_tot);
      fprintf(file1, " %i %i %6.3e %6.3e \n", nside, L, accuracy_tot,timing_tot);

      accuracy_tot = 0.0;
      timing_tot = 0.0;
      for (repeat=0; repeat<NREPEAT; repeat++){
        s2let_axisym_hpx_wav_test(&accuracy, &timing, nside, B, L, J_min, seed);
        accuracy_tot += accuracy;
        timing_tot += timing;
      }
      accuracy_tot /= NREPEAT;
      timing_tot /= NREPEAT;
      printf("2-WAV: Nside = %i ; Lmax = %i : Accuracy : %6.3e in %6.3e sec \n",nside, L, accuracy_tot,timing_tot);
      fprintf(file2, " %i %i %6.3e %6.3e \n", nside, L, accuracy_tot,timing_tot);

    }
  } 

  close(file1);
  close(file2);
*/
  


  printf("==========================================================\n");
  return 0;		
}

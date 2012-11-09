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

/*!
 * PROGRAM : s2let_hpx_demo
 * COMMAND : bin/s2let_hpx_demo
 * ARGUMENTS : none
 */
int main(int argc, char *argv[]) 
{
  printf("--------------------------------------------------\n");
  printf(" S2LET : Healpix wavelet transform \n");
  printf("--------------------------------------------------\n");

  const int L = 256;    // Harmonic band-limit
  const int B = 4;      // Wavelet parameters
  const int J_min = 2;  // First wavelet scale to use

  // The input file is a random CMB simulation (healpix map with nside=128)
  char file[100] = "data/somecmbsimu_hpx_128.fits";
  const int nside = s2let_fits_hpx_read_nside(file);

  // Read Healpix map from file
  double *f = (double*)calloc(12*nside*nside, sizeof(double));
  s2let_hpx_read_map(f, file, nside);

  // Allocate space for wavelet maps (corresponding to the triplet B/L/J_min)
  double *f_wav, *f_scal;
  s2let_axisym_hpx_allocate_f_wav_real(&f_wav, &f_scal, nside, B, L, J_min);

  // Perform wavelet analysis from scratch with all signals given as Healpix maps
  clock_t time_start = clock();
  s2let_axisym_hpx_wav_analysis_real(f_wav, f_scal, f, nside, B, L, J_min);
  clock_t time_end = clock();
  printf(" - Wavelet analysis   : %4.4f seconds\n", 
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);
	
  // Output the wavelets to FITS files
  char outfile[100];
  char params[100];
  sprintf(params, "%d%s%d%s%d", L, "_", B, "_", J_min);
  int j, J = s2let_j_max(L, B); // Explicitly compute the maximum wavelet scale
  int offset = 0; // Start with the first wavelet
  for(j = J_min; j <= J; j++){
    sprintf(outfile, "%s%s%s%s%d%s", "data/somecmbsimu_hpx_128", "_wav_", params, "_", j, ".fits");
    printf(" Outfile_wav[j=%i] = %s\n",j,outfile);
    remove(outfile); // In case the file exists
    s2let_hpx_write_map(outfile, f_wav + offset, nside); // Now write the map to fits file
    offset += 12 * nside * nside; // Go to the next wavelet
  }
  // Finally write the scaling function
  sprintf(outfile, "%s%s%s%s", "data/somecmbsimu_hpx_128", "_scal_", params, ".fits");
  printf(" Outfile_scal = %s\n",outfile);
  remove(outfile); // In case the file exists
  s2let_hpx_write_map(outfile, f_scal, nside); // Now write the map to fits file

  free(f_wav);
  free(f_scal);
  free(f);

  printf("--------------------------------------------------\n");

  return 0;		
}


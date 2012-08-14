// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <assert.h>

int main(int argc, char *argv[]) 
{
	
	const int L = 128;    // Harmonic band-limit
	const int B = 3;      // Wavelet parameters
	const int J_min = 1;  // First wavelet scale to use

	// The input file is a random CMB simulation (healpix map with nside=128)
	char file[100] = "data/some_cmb_simu.fits";
	const int nside = 128;

	// Read Healpix map from file
	double *f = (double*)calloc(12*nside*nside, sizeof(double));
	read_healpix_map(f, file, nside);

	// Allocate space for wavelet maps (corresponding to the triplet B/L/J_min)
	double *f_wav, *f_scal;
	s2let_axisym_hpx_allocate_f_wav_real(&f_wav, &f_scal, nside, B, L, J_min);

	// Perform wavelet analysis from scratch with all signals given as Healpix maps
	clock_t time_start = clock();
	s2let_axisym_hpx_wav_analysis_real(f_wav, f_scal, f, nside, B, L, J_min);
	clock_t time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	
	// Output the wavelets to FITS files
    char outfile[100];
    char params[100];
    sprintf(params, "%d%s%d%s%d", L, "_", B, "_", J_min);
	int j, J = s2let_j_max(L, B); // Explicitly compute the maximum wavelet scale
	int offset = 0; // Start with the first wavelet
	for(j = J_min; j <= J; j++){
		sprintf(outfile, "%s%s%s%s%d%s", "some_cmb_simu", "_wav_", params, "_", j, ".fits");
		printf("Outfile_wav[j=%i] = %s\n",j,outfile);
		remove(outfile); // In case the file exists
		write_healpix_map(outfile, f_wav + offset, nside); // Now write the map to fits file
		offset += 12 * nside * nside; // Go to the next wavelet
	}
	// Finally write the scaling function
	sprintf(outfile, "%s%s%s%s", "some_cmb_simu", "_scal_", params, ".fits");
	printf("Outfile_scal = %s\n",outfile);
	remove(outfile); // In case the file exists
	write_healpix_map(outfile, f_scal, nside); // Now write the map to fits file

	return 0;		
}



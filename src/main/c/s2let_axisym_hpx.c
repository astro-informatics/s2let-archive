// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <complex.h>
#include <stdlib.h>

void s2let_axisym_hpx_allocate_f_wav_real(double **f_wav, double **f_scal, int nside, int B, int L, int J_min)
{
	int J = s2let_j_max(L, B);
	*f_wav = (double*)calloc((J+1-J_min) * 12*nside*nside, sizeof(double));
	*f_scal = (double*)calloc(12*nside*nside, sizeof(double));
}


void s2let_axisym_hpx_wav_analysis_real(double *f_wav, double *f_scal, const double *f, int nside, int B, int L, int J_min)
{
	int bandlimit, j, offset, offset_lm;
	int J = s2let_j_max(L, B);

	double *wav_lm, *scal_lm;
	s2let_axisym_allocate_wav_lm(&wav_lm, &scal_lm, B, L);
	s2let_axisym_wav_lm(wav_lm, scal_lm, B, L, J_min);

	complex double *flm, *f_wav_lm, *f_scal_lm;
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	s2let_axisym_allocate_f_wav_multires_lm(&f_wav_lm, &f_scal_lm, B, L, J_min);

	s2let_hpx_map2alm_real(flm, f, nside, L);

	s2let_axisym_wav_analysis_multires_lm(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, B, L, J_min);

	bandlimit = MIN(s2let_bandlimit(B, J_min-1), L);
	s2let_hpx_alm2map_real(f_scal, f_scal_lm, nside, bandlimit);

	offset = 0;
	offset_lm = 0;
	for(j = J_min; j <= J; j++){
		bandlimit = MIN(s2let_bandlimit(B, j), L);
		s2let_hpx_alm2map_real(f_wav + offset, f_wav_lm + offset_lm, nside, bandlimit);
		offset_lm += bandlimit * bandlimit;
		offset += 12 * nside * nside;
	}

	free(flm);
	free(f_scal_lm);
	free(f_wav_lm);
}

void s2let_axisym_hpx_wav_synthesis_real(double *f, const double *f_wav, const double *f_scal, int nside, int B, int L, int J_min)
{
	int bandlimit, j, offset, offset_lm;
	int J = s2let_j_max(L, B);

	double *wav_lm, *scal_lm;
	s2let_axisym_allocate_wav_lm(&wav_lm, &scal_lm, B, L);
	s2let_axisym_wav_lm(wav_lm, scal_lm, B, L, J_min);

	complex double *flm, *f_wav_lm, *f_scal_lm;
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	s2let_axisym_allocate_f_wav_multires_lm(&f_wav_lm, &f_scal_lm, B, L, J_min);

	bandlimit = MIN(s2let_bandlimit(B, J_min-1), L);
	s2let_hpx_map2alm_real(f_scal_lm, f_scal, nside, bandlimit);
	offset = 0;
	offset_lm = 0;
	for(j = J_min; j <= J; j++){
		bandlimit = MIN(s2let_bandlimit(B, j), L);
		s2let_hpx_map2alm_real(f_wav_lm + offset_lm, f_wav + offset, nside, bandlimit);
		offset_lm += bandlimit * bandlimit;
		offset += 12 * nside * nside;
	}

	s2let_axisym_wav_synthesis_multires_lm(flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, B, L, J_min);

	s2let_hpx_alm2map_real(f, flm, nside, L);

	free(flm);
	free(f_scal_lm);
	free(f_wav_lm);
}

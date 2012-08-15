// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"

/*!
 * Allocates arrays for wavelets and scaling functions in pixel space (Healpix sampling).
 *
 * \param[out]  f_wav Array of wavelets HEALPIX maps.
 * \param[out]  f_scal Scaling function HEALPIX map.
 * \param[in]  nside HEALPIX resolution parameter.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_hpx_allocate_f_wav_real(double **f_wav, double **f_scal, int nside, int B, int L, int J_min)
{
	int J = s2let_j_max(L, B);
	*f_wav = (double*)calloc((J+1-J_min) * 12*nside*nside, sizeof(double));
	*f_scal = (double*)calloc(12*nside*nside, sizeof(double));
}


/*!
 * Perform wavelet transform in real space (from scratch, gives wavelet maps).
 * Input function is real.
 * Sampling scheme : HEALPIX sampling.
 * Spherical wavelets : analysis in real space, HEALPIX sampling.
 * Note : multiresolution in used in harmonic space but all maps are at resolution nside.
 *
 * \param[out]  f_wav Array of wavelets HEALPIX maps.
 * \param[out]  f_scal Scaling function HEALPIX map.
 * \param[in]  f Input function (HEALPIX map)
 * \param[in]  nside HEALPIX resolution parameter.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
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

	healpix_forward_real(flm, f, nside, L);

	s2let_axisym_wav_analysis_multires_lm(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, B, L, J_min);

	bandlimit = MIN(s2let_bandlimit(B, J_min-1), L);
	healpix_inverse_real(f_scal, f_scal_lm, nside, bandlimit);

	offset = 0;
	offset_lm = 0;
	for(j = J_min; j <= J; j++){
		bandlimit = MIN(s2let_bandlimit(B, j), L);
		healpix_inverse_real(f_wav + offset, f_wav_lm + offset_lm, nside, bandlimit);
		offset_lm += bandlimit * bandlimit;
		offset += 12 * nside * nside;
	}
	offset_lm -= bandlimit * bandlimit;
	offset -= 12 * nside * nside;

	free(flm);
	free(f_scal_lm);
	free(f_wav_lm);
}

/*!
 * Perform wavelet transform in real space (from scratch, gives wavelet maps).
 * Input function is real.
 * Sampling scheme : HEALPIX sampling.
 * Spherical wavelets : synthesis in real space, HEALPIX sampling.
 * Note : multiresolution in used in harmonic space but all maps are at resolution nside.
 *
 * \param[out]  f Input function (HEALPIX map)
 * \param[in]  f_wav Array of wavelets HEALPIX maps.
 * \param[in]  f_scal Scaling function HEALPIX map.
 * \param[in]  nside HEALPIX resolution parameter.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
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
	healpix_forward_real(f_scal_lm, f_scal, nside, bandlimit);
	offset = 0;
	offset_lm = 0;
	for(j = J_min; j <= J; j++){
		bandlimit = MIN(s2let_bandlimit(B, j), L);
		healpix_forward_real(f_wav_lm + offset_lm, f_wav + offset, nside, bandlimit);
		offset_lm += bandlimit * bandlimit;
		offset += 12 * nside * nside;
	}

	s2let_axisym_wav_synthesis_multires_lm(flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, B, L, J_min);

	healpix_inverse_real(f, flm, nside, L);

	free(flm);
	free(f_scal_lm);
	free(f_wav_lm);
}

// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <complex.h> 
#include <ssht.h>
#include <stdlib.h>

/*!
 * Allocates arrays for wavelets and scaling functions in pixel space (MW sampling).
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_allocate_f_wav_multires(complex double **f_wav, complex double **f_scal, int B, int L, int J_min)
{
	int J = s2let_j_max(L, B);
	int j, bandlimit, total = 0;
	for(j = J_min; j <= J; j++){
		bandlimit = MIN(s2let_bandlimit(B, j), L);
		total += bandlimit * (2 * bandlimit - 1);
	}
	*f_wav = (complex double*)calloc(total, sizeof(complex double));
	bandlimit = MIN(s2let_bandlimit(B, J_min-1), L);
	*f_scal = (complex double*)calloc(bandlimit * (2*bandlimit-1), sizeof(complex double));
}

/*!
 * Allocates arrays for multiresolution wavelets and scaling functions in pixel space (MW sampling).
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_allocate_f_wav_multires_real(double **f_wav, double **f_scal, int B, int L, int J_min)
{
	int J = s2let_j_max(L, B);
	int j, bandlimit, total = 0;
	for(j = J_min; j <= J; j++){
		bandlimit = MIN(s2let_bandlimit(B, j), L);
		total += bandlimit * (2 * bandlimit - 1);
	}
	*f_wav = (double*)calloc(total, sizeof(double));
	bandlimit = MIN(s2let_bandlimit(B, J_min-1), L);
	*f_scal = (double*)calloc(bandlimit * (2*bandlimit-1), sizeof(double));
}

/*!
 * Allocates arrays for final wavelets and scaling functions in pixel space (MW sampling).
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_allocate_f_wav(complex double **f_wav, complex double **f_scal, int B, int L, int J_min)
{
	int J = s2let_j_max(L, B);
	*f_wav = (complex double*)calloc((J+1-J_min) * L *(2*L-1), sizeof(complex double));
	*f_scal = (complex double*)calloc(L * (2*L-1), sizeof(complex double));
}

/*!
 * Allocates arrays for final wavelets and scaling functions in pixel space (MW sampling).
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_allocate_f_wav_real(double **f_wav, double **f_scal, int B, int L, int J_min)
{
	int J = s2let_j_max(L, B);
	*f_wav = (double*)calloc((J+1-J_min) * L *(2*L-1), sizeof(double));
	*f_scal = (double*)calloc(L * (2*L-1), sizeof(double));
}

/*!
 * Spherical wavelets : full resolution analysis in real space, MW sampling.
 * Perform wavelet transform in real space (from scratch, gives wavelet maps).
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling..
 * \param[in]  f Input function (MW sampling)
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_wav_analysis(complex double *f_wav, complex double *f_scal, const complex double *f, int B, int L, int J_min)
{
	int spin = 0;
	int verbosity = 0;
	ssht_dl_method_t dl_method = SSHT_DL_RISBO;

	int j, offset, offset_lm;
	int J = s2let_j_max(L, B);
	//int l_min = s2let_axisym_el_min(B, J_min);

	double *wav_lm, *scal_lm;
	s2let_axisym_allocate_wav_lm(&wav_lm, &scal_lm, B, L);
	s2let_axisym_wav_lm(wav_lm, scal_lm, B, L, J_min);

	complex double *flm, *f_wav_lm, *f_scal_lm;
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	s2let_axisym_allocate_f_wav_lm(&f_wav_lm, &f_scal_lm, B, L, J_min);

	ssht_core_mw_forward_sov_conv_sym(flm, f, L, spin, dl_method, verbosity);

	s2let_axisym_wav_analysis_lm(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, B, L, J_min);

	ssht_core_mw_inverse_sov_sym(f_scal, f_scal_lm, L, spin, dl_method, verbosity);
	offset = 0;
	offset_lm = 0;
	for(j = J_min; j <= J; j++){
		ssht_core_mw_inverse_sov_sym(f_wav + offset, f_wav_lm + offset_lm, L, spin, dl_method, verbosity);
		offset_lm += L * L;
		offset += L * (2 * L - 1);
	}

	free(flm);
	free(f_scal_lm);
	free(f_wav_lm);
}

/*!
 * Perform wavelet transform in real space (from scratch, gives wavelet maps).
 * Sampling scheme : MW sampling.
 * Spherical wavelets : synthesis in real space, MW sampling.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Array of wavelets maps, MW sampling.
 * \param[in]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_wav_synthesis(complex double *f, const complex double *f_wav, const complex double *f_scal, int B, int L, int J_min)
{
	int spin = 0;
	int verbosity = 0;
	ssht_dl_method_t dl_method = SSHT_DL_RISBO;

	int j, offset, offset_lm;
	int J = s2let_j_max(L, B);
	//int l_min = s2let_axisym_el_min(B, J_min);

	double *wav_lm, *scal_lm;
	s2let_axisym_allocate_wav_lm(&wav_lm, &scal_lm, B, L);
	s2let_axisym_wav_lm(wav_lm, scal_lm, B, L, J_min);

	complex double *flm, *f_wav_lm, *f_scal_lm;
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	s2let_axisym_allocate_f_wav_lm(&f_wav_lm, &f_scal_lm, B, L, J_min);

	ssht_core_mw_forward_sov_conv_sym(f_scal_lm, f_scal, L, spin, dl_method, verbosity);
	offset = 0;
	offset_lm = 0;
	for(j = J_min; j <= J; j++){
		ssht_core_mw_forward_sov_conv_sym(f_wav_lm + offset_lm, f_wav + offset, L, spin, dl_method, verbosity);
		offset_lm += L * L;
		offset += L * (2 * L - 1);
	}

	s2let_axisym_wav_synthesis_lm(flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, B, L, J_min);

	ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

	free(flm);
	free(f_scal_lm);
	free(f_wav_lm);
}


/*!
 * Perform wavelet transform in real space (from scratch, gives wavelet maps).
 * Input function is real.
 * Sampling scheme : MW sampling.
 * Spherical wavelets : analysis in real space, MW sampling.
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  f Input function (MW sampling)
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_wav_analysis_real(double *f_wav, double *f_scal, const double *f, int B, int L, int J_min)
{
	int verbosity = 0;
	ssht_dl_method_t dl_method = SSHT_DL_RISBO;

	int j, offset, offset_lm;
	int J = s2let_j_max(L, B);
	//int l_min = s2let_axisym_el_min(B, J_min);

	double *wav_lm, *scal_lm;
	s2let_axisym_allocate_wav_lm(&wav_lm, &scal_lm, B, L);
	s2let_axisym_wav_lm(wav_lm, scal_lm, B, L, J_min);

	complex double *flm, *f_wav_lm, *f_scal_lm;
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	s2let_axisym_allocate_f_wav_lm(&f_wav_lm, &f_scal_lm, B, L, J_min);

	ssht_core_mw_forward_sov_conv_sym_real(flm, f, L, dl_method, verbosity);

	s2let_axisym_wav_analysis_lm(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, B, L, J_min);

	ssht_core_mw_inverse_sov_sym_real(f_scal, f_scal_lm, L, dl_method, verbosity);
	offset = 0;
	offset_lm = 0;
	for(j = J_min; j <= J; j++){
		ssht_core_mw_inverse_sov_sym_real(f_wav + offset, f_wav_lm + offset_lm, L, dl_method, verbosity);
		offset_lm += L * L;
		offset += L * (2 * L - 1);
	}

	free(flm);
	free(f_scal_lm);
	free(f_wav_lm);
}

/*!
 * Perform wavelet transform in real space (from scratch, gives wavelet maps).
 * Input function is real.
 * Sampling scheme : MW sampling.
 * Spherical wavelets : synthesis in real space, MW sampling.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Array of wavelets maps, MW sampling.
 * \param[in]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_wav_synthesis_real(double *f, const double *f_wav, const double *f_scal, int B, int L, int J_min)
{
	int verbosity = 0;
	ssht_dl_method_t dl_method = SSHT_DL_RISBO;

	int j, offset, offset_lm;
	int J = s2let_j_max(L, B);
	//int l_min = s2let_axisym_el_min(B, J_min);

	double *wav_lm, *scal_lm;
	s2let_axisym_allocate_wav_lm(&wav_lm, &scal_lm, B, L);
	s2let_axisym_wav_lm(wav_lm, scal_lm, B, L, J_min);

	complex double *flm, *f_wav_lm, *f_scal_lm;
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	s2let_axisym_allocate_f_wav_lm(&f_wav_lm, &f_scal_lm, B, L, J_min);

	ssht_core_mw_forward_sov_conv_sym_real(f_scal_lm, f_scal, L, dl_method, verbosity);
	offset = 0;
	offset_lm = 0;
	for(j = J_min; j <= J; j++){
		ssht_core_mw_forward_sov_conv_sym_real(f_wav_lm + offset_lm, f_wav + offset, L, dl_method, verbosity);
		offset_lm += L * L;
		offset += L * (2 * L - 1);
	}

	s2let_axisym_wav_synthesis_lm(flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, B, L, J_min);

	ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

	free(flm);
	free(f_scal_lm);
	free(f_wav_lm);
}


/*!
 * Spherical wavelets : multiresolution analysis in real space, MW sampling.
 * Perform multiresolution wavelet transform in real space (from scratch, gives wavelet maps).
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  f Input function (MW sampling)
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_wav_analysis_multires(complex double *f_wav, complex double *f_scal, const complex double *f, int B, int L, int J_min)
{
	int spin = 0;
	int verbosity = 0;
	ssht_dl_method_t dl_method = SSHT_DL_RISBO;

	int bandlimit, j, offset, offset_lm;
	int J = s2let_j_max(L, B);
	//int l_min = s2let_axisym_el_min(B, J_min);

	double *wav_lm, *scal_lm;
	s2let_axisym_allocate_wav_lm(&wav_lm, &scal_lm, B, L);
	s2let_axisym_wav_lm(wav_lm, scal_lm, B, L, J_min);

	complex double *flm, *f_wav_lm, *f_scal_lm;
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	s2let_axisym_allocate_f_wav_multires_lm(&f_wav_lm, &f_scal_lm, B, L, J_min);

	ssht_core_mw_forward_sov_conv_sym(flm, f, L, spin, dl_method, verbosity);

	s2let_axisym_wav_analysis_multires_lm(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, B, L, J_min);

	bandlimit = MIN(s2let_bandlimit(B, J_min-1), L);
	ssht_core_mw_inverse_sov_sym(f_scal, f_scal_lm, bandlimit, spin, dl_method, verbosity);
	offset = 0;
	offset_lm = 0;
	for(j = J_min; j <= J; j++){
		bandlimit = MIN(s2let_bandlimit(B, j), L);
		ssht_core_mw_inverse_sov_sym(f_wav + offset, f_wav_lm + offset_lm, bandlimit, spin, dl_method, verbosity);
		offset_lm += bandlimit * bandlimit;
		offset += bandlimit * (2 * bandlimit - 1);
	}

	free(flm);
	free(f_scal_lm);
	free(f_wav_lm);
}


/*!
 * Spherical wavelets : multiresolution synthesis in real space, MW sampling.
 * Perform multiresolution wavelet transform in real space (from scratch, gives wavelet maps).
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Array of wavelets maps, MW sampling.
 * \param[in]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_wav_synthesis_multires(complex double *f, const complex double *f_wav, const complex double *f_scal, int B, int L, int J_min)
{
	int spin = 0;
	int verbosity = 0;
	ssht_dl_method_t dl_method = SSHT_DL_RISBO;

	int bandlimit, j, offset, offset_lm;
	int J = s2let_j_max(L, B);
	//int l_min = s2let_axisym_el_min(B, J_min);

	double *wav_lm, *scal_lm;
	s2let_axisym_allocate_wav_lm(&wav_lm, &scal_lm, B, L);
	s2let_axisym_wav_lm(wav_lm, scal_lm, B, L, J_min);

	complex double *flm, *f_wav_lm, *f_scal_lm;
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	s2let_axisym_allocate_f_wav_multires_lm(&f_wav_lm, &f_scal_lm, B, L, J_min);

	bandlimit = MIN(s2let_bandlimit(B, J_min-1), L);
	ssht_core_mw_forward_sov_conv_sym(f_scal_lm, f_scal, bandlimit, spin, dl_method, verbosity);
	offset = 0;
	offset_lm = 0;
	for(j = J_min; j <= J; j++){
		bandlimit = MIN(s2let_bandlimit(B, j), L);
		ssht_core_mw_forward_sov_conv_sym(f_wav_lm + offset_lm, f_wav + offset, bandlimit, spin, dl_method, verbosity);
		offset_lm += bandlimit * bandlimit;
		offset += bandlimit * (2 * bandlimit - 1);
	}

	s2let_axisym_wav_synthesis_multires_lm(flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, B, L, J_min);

	ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

	free(flm);
	free(f_scal_lm);
	free(f_wav_lm);
}


/*!
 * Spherical wavelets : multiresolution analysis in real space, MW sampling.
 * Perform multiresolution wavelet transform in real space (from scratch, gives wavelet maps).
 * Input function is real.
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  f Input function (MW sampling)
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_wav_analysis_multires_real(double *f_wav, double *f_scal, const double *f, int B, int L, int J_min)
{
	int verbosity = 0;
	ssht_dl_method_t dl_method = SSHT_DL_RISBO;

	int bandlimit, j, offset, offset_lm;
	int J = s2let_j_max(L, B);
	//int l_min = s2let_axisym_el_min(B, J_min);

	double *wav_lm, *scal_lm;
	s2let_axisym_allocate_wav_lm(&wav_lm, &scal_lm, B, L);
	s2let_axisym_wav_lm(wav_lm, scal_lm, B, L, J_min);

	complex double *flm, *f_wav_lm, *f_scal_lm;
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	s2let_axisym_allocate_f_wav_multires_lm(&f_wav_lm, &f_scal_lm, B, L, J_min);

	ssht_core_mw_forward_sov_conv_sym_real(flm, f, L, dl_method, verbosity);

	s2let_axisym_wav_analysis_multires_lm(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, B, L, J_min);

	bandlimit = MIN(s2let_bandlimit(B, J_min-1), L);
	ssht_core_mw_inverse_sov_sym_real(f_scal, f_scal_lm, bandlimit, dl_method, verbosity);
	offset = 0;
	offset_lm = 0;
	for(j = J_min; j <= J; j++){
		bandlimit = MIN(s2let_bandlimit(B, j), L);
		ssht_core_mw_inverse_sov_sym_real(f_wav + offset, f_wav_lm + offset_lm, bandlimit, dl_method, verbosity);
		offset_lm += bandlimit * bandlimit;
		offset += bandlimit * (2 * bandlimit - 1);
	}

	free(flm);
	free(f_scal_lm);
	free(f_wav_lm);
}

/*!
 * Spherical wavelets : multiresolution synthesis in real space, MW sampling.
 * Perform wavelet transform in real space (from scratch, gives wavelet maps).
 * Input function is real.
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Array of wavelets maps, MW sampling.
 * \param[in]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_wav_synthesis_multires_real(double *f, const double *f_wav, const double *f_scal, int B, int L, int J_min)
{
	int verbosity = 0;
	ssht_dl_method_t dl_method = SSHT_DL_RISBO;

	int bandlimit, j, offset, offset_lm;
	int J = s2let_j_max(L, B);
	//int l_min = s2let_axisym_el_min(B, J_min);

	double *wav_lm, *scal_lm;
	s2let_axisym_allocate_wav_lm(&wav_lm, &scal_lm, B, L);
	s2let_axisym_wav_lm(wav_lm, scal_lm, B, L, J_min);

	complex double *flm, *f_wav_lm, *f_scal_lm;
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	s2let_axisym_allocate_f_wav_multires_lm(&f_wav_lm, &f_scal_lm, B, L, J_min);

	bandlimit = MIN(s2let_bandlimit(B, J_min-1), L);
	ssht_core_mw_forward_sov_conv_sym_real(f_scal_lm, f_scal, bandlimit, dl_method, verbosity);
	offset = 0;
	offset_lm = 0;
	for(j = J_min; j <= J; j++){
		bandlimit = MIN(s2let_bandlimit(B, j), L);
		ssht_core_mw_forward_sov_conv_sym_real(f_wav_lm + offset_lm, f_wav + offset, bandlimit, dl_method, verbosity);
		offset_lm += bandlimit * bandlimit;
		offset += bandlimit * (2 * bandlimit - 1);
	}

	s2let_axisym_wav_synthesis_multires_lm(flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, B, L, J_min);

	ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

	free(flm);
	free(f_scal_lm);
	free(f_wav_lm);
}


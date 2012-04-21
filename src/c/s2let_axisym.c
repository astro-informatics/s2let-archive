// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))

/*!
 * Allocates arrays for multiresolution wavelets and scaling functions in pixel space.
 *
 * \param[out]  f_wav Wavelets.
 * \param[out]  f_scal Scaling function.
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
 * Allocates arrays for multiresolution wavelets and scaling functions in pixel space.
 *
 * \param[out]  f_wav Wavelets.
 * \param[out]  f_scal Scaling function.
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
 * Allocates arrays for final wavelets and scaling functions in pixel space.
 *
 * \param[out]  f_wav Wavelets.
 * \param[out]  f_scal Scaling function.
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
 * Allocates arrays for final wavelets and scaling functions in pixel space.
 *
 * \param[out]  f_wav Wavelets.
 * \param[out]  f_scal Scaling function.
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
 * Allocates arrays for final wavelets and scaling functions in harmonic space.
 *
 * \param[out]  f_wav_lm Wavelets.
 * \param[out]  f_scal_lm Scaling function.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_axisym_allocate_f_wav_lm(complex double **f_wav_lm, complex double **f_scal_lm, int B, int L, int J_min)
{
	int J = s2let_j_max(L, B);
	*f_wav_lm = (complex double*)calloc((J+1-J_min) * L * L, sizeof(complex double));
	*f_scal_lm = (complex double*)calloc(L * L, sizeof(complex double));
}

/*!
 * Allocates multiresolution arrays for final wavelets and scaling functions in harmonic space.
 *
 * \param[out]  f_wav_lm Wavelets.
 * \param[out]  f_scal_lm Scaling function.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_axisym_allocate_f_wav_multires_lm(complex double **f_wav_lm, complex double **f_scal_lm, int B, int L, int J_min)
{
	int J = s2let_j_max(L, B);
	int j, bandlimit, total = 0;
	for(j = J_min; j <= J; j++){
		bandlimit = MIN(s2let_bandlimit(B, j), L);
		total += bandlimit * bandlimit;
	}
	*f_wav_lm = (complex double*)calloc(total, sizeof(complex double));
	bandlimit = MIN(s2let_bandlimit(B, J_min-1), L);
	*f_scal_lm = (complex double*)calloc(bandlimit * bandlimit, sizeof(complex double));
}


/*!
 * Allocates arrays for the kernels of the wavelets and the scaling functions.
 *
 * \param[out]  wav_lm Wavelet kernels.
 * \param[out]  scal_lm Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_axisym_allocate_wav_lm(double **wav_lm, double **scal_lm, int B, int L)
{
	int J = s2let_j_max(L, B);
	*wav_lm = (double*)calloc((J+1) * L, sizeof(double));
	*scal_lm = (double*)calloc(L, sizeof(double));
}

/*!
 * Computes the kernels of the wavelets and the scaling functions.
 *
 * \param[out]  wav_lm Wavelet kernels.
 * \param[out]  scal_lm Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_wav_lm(double *wav_lm, double *scal_lm, int B, int L, int J_min)
{
	int j, l;
	int J = s2let_j_max(L, B);
	//int J_min = 0;
	//int l_min = s2let_axisym_el_min(B, J_min);
	double k0;
	double *kappa, *kappa0;
	s2let_axisym_allocate_tilling(&kappa, &kappa0, B, L);
	s2let_axisym_tilling(kappa, kappa0, B, L, J_min);

	for (j = J_min; j <= J; j++){
		for (l = 0; l < L; l++){
			k0 = sqrt( (2 * l + 1) / (4.0 * PI) ) * kappa[l+j*L];
			wav_lm[j*L+l] = k0;
		}
	}
	for (l = 0; l < L; l++){
		k0 = sqrt( (2 * l + 1) / (4.0 * PI) ) * kappa0[l];
		scal_lm[l] = k0;
	}

	free(kappa);
	free(kappa0);
}

/*!
 * Perform wavelet transform in harmonic space (from precomputed kernels, gives SHA coefficients).
 * Spherical wavelets : analysis in harmonic space.
 *
 * \param[out]  f_wav_lm Wavelet transform (SHA of wavelet contribution).
 * \param[out]  f_scal_lm Wavelet transform (SHA of scaling contribution).
 * \param[in]  flm Spherical harmonic decomposition of input function.
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_lm Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_wav_analysis_lm(complex double *f_wav_lm, complex double *f_scal_lm, const complex double *flm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min)
{
	int offset, j, l, m;
	int J = s2let_j_max(L, B);
	double wav0, scal0;
	//int l_min = s2let_axisym_el_min(B, J_min);

	offset = 0;
	for (j = J_min; j <= J; j++){
		for (l = 0; l < L; l++){
			wav0 = sqrt((4.0*PI)/(2.0*l+1.0)) * wav_lm[j*L+l];
			for (m = -l; m <= l; m++){
				f_wav_lm[offset + l*l + l + m] = flm[lm2ind(l,m)] * wav0 ;
			}
		}
		offset += L * L;
	}
	for (l = 0; l < L; l++){
		scal0 = sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lm[l];
		for (m = -l; m <= l; m++){
			f_scal_lm[lm2ind(l,m)] = flm[lm2ind(l,m)] * scal0 ;
		}
	}
}

/*!
 * Perform wavelet transform in harmonic space (from precomputed kernels, gives SHA coefficients).
 * Spherical wavelets : synthesis in harmonic space.
 *
 * \param[out]  flm Spherical harmonic decomposition of input function.
 * \param[in]  f_wav_lm Wavelet transform (SHA of wavelet contribution).
 * \param[in]  f_scal_lm Wavelet transform (SHA of scaling contribution).
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_lm Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_wav_synthesis_lm(complex double *flm, const complex double *f_wav_lm, const complex double *f_scal_lm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min)
{
	int offset, j, l, m;
	int J = s2let_j_max(L, B);
	double wav0, scal0;
	//int l_min = s2let_axisym_el_min(B, J_min);

	offset = 0; 
	for (j = J_min; j <= J; j++){
		for (l = 0; l < L; l++){
			wav0 = sqrt((4.0*PI)/(2.0*l+1.0)) * wav_lm[j*L+l];
			for (m = -l; m <= l; m++){
				flm[lm2ind(l,m)] += f_wav_lm[offset + l*l + l + m] * wav0 ;
			}
		}
		offset += L * L;
	}
	for (l = 0; l < L; l++){
		scal0 = sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lm[l];
		for (m = -l; m <= l; m++){
			flm[lm2ind(l,m)] += f_scal_lm[lm2ind(l,m)] * scal0 ;
		}
	}
}

/*!
 * Perform multiresolution wavelet transform in harmonic space (from precomputed kernels, gives SHA coefficients).
 * Spherical wavelets : analysis in harmonic space.
 *
 * \param[out]  f_wav_lm Wavelet transform (SHA of wavelet contribution).
 * \param[out]  f_scal_lm Wavelet transform (SHA of scaling contribution).
 * \param[in]  flm Spherical harmonic decomposition of input function.
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_lm Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_wav_analysis_multires_lm(complex double *f_wav_lm, complex double *f_scal_lm, const complex double *flm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min)
{
	int bandlimit, offset, j, l, m;
	int J = s2let_j_max(L, B);
	double wav0, scal0;
	//int l_min = s2let_axisym_el_min(B, J_min);

	offset = 0;
	for (j = J_min; j <= J; j++){
		bandlimit = MIN(s2let_bandlimit(B, j), L);
		for (l = 0; l < bandlimit; l++){
			wav0 = sqrt((4.0*PI)/(2.0*l+1.0)) * wav_lm[j*L+l];
			for (m = -l; m <= l; m++){
				f_wav_lm[offset + l*l + l + m] = flm[lm2ind(l,m)] * wav0 ;
			}
		}
		offset += bandlimit * bandlimit;
	}
	bandlimit = MIN(s2let_bandlimit(B, J_min-1), L);
	for (l = 0; l < bandlimit; l++){
		scal0 = sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lm[l];
		for (m = -l; m <= l; m++){
			f_scal_lm[lm2ind(l,m)] = flm[lm2ind(l,m)] * scal0 ;
		}
	}
}

/*!
 * Perform multiresolution wavelet transform in harmonic space (from precomputed kernels, gives SHA coefficients).
 * Spherical wavelets : synthesis in harmonic space.
 *
 * \param[out]  flm Spherical harmonic decomposition of input function.
 * \param[in]  f_wav_lm Wavelet transform (SHA of wavelet contribution).
 * \param[in]  f_scal_lm Wavelet transform (SHA of scaling contribution).
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_lm Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_wav_synthesis_multires_lm(complex double *flm, const complex double *f_wav_lm, const complex double *f_scal_lm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min)
{
	int bandlimit, offset, j, l, m;
	int J = s2let_j_max(L, B);
	double wav0, scal0;
	//int l_min = s2let_axisym_el_min(B, J_min);

	offset = 0; 
	for (j = J_min; j <= J; j++){
		bandlimit = MIN(s2let_bandlimit(B, j), L);
		for (l = 0; l < bandlimit; l++){
			wav0 = sqrt((4.0*PI)/(2.0*l+1.0)) * wav_lm[j*L+l];
			for (m = -l; m <= l; m++){
				flm[lm2ind(l,m)] += f_wav_lm[offset + l*l + l + m] * wav0 ;
			}
		}
		offset += bandlimit * bandlimit;
	}
	bandlimit = MIN(s2let_bandlimit(B, J_min-1), L);
	for (l = 0; l < bandlimit; l++){
		scal0 = sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lm[l];
		for (m = -l; m <= l; m++){
			flm[lm2ind(l,m)] += f_scal_lm[lm2ind(l,m)] * scal0 ;
		}
	}
}

/*!
 * Perform wavelet transform in real space (from scratch, gives pixel space components).
 * Sampling scheme : MW sampling.
 * Spherical wavelets : analysis in real space.
 *
 * \param[out]  f_wav Wavelet transform (wavelet contribution in real space).
 * \param[out]  f_scal Wavelet transform (scaling contribution in real space).
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
 * Perform wavelet transform in real space (from scratch, gives pixel space components).
 * Sampling scheme : MW sampling.
 * Spherical wavelets : synthesis in real space.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Wavelet transform (wavelet contribution in real space).
 * \param[in]  f_scal Wavelet transform (scaling contribution in real space).
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
 * Perform wavelet transform in real space (from scratch, gives pixel space components).
 * Input function is real.
 * Sampling scheme : MW sampling.
 * Spherical wavelets : analysis in real space.
 *
 * \param[out]  f_wav Wavelet transform (wavelet contribution in real space).
 * \param[out]  f_scal Wavelet transform (scaling contribution in real space).
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
 * Perform wavelet transform in real space (from scratch, gives pixel space components).
 * Input function is real.
 * Sampling scheme : MW sampling.
 * Spherical wavelets : synthesis in real space.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Wavelet transform (wavelet contribution in real space).
 * \param[in]  f_scal Wavelet transform (scaling contribution in real space).
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
 * Perform multiresolution wavelet transform in real space (from scratch, gives pixel space components).
 * Sampling scheme : MW sampling.
 * Spherical wavelets : analysis in real space.
 *
 * \param[out]  f_wav Wavelet transform (wavelet contribution in real space).
 * \param[out]  f_scal Wavelet transform (scaling contribution in real space).
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
 * Perform wavelet transform in real space (from scratch, gives pixel space components).
 * Sampling scheme : MW sampling.
 * Spherical wavelets : synthesis in real space.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Wavelet transform (wavelet contribution in real space).
 * \param[in]  f_scal Wavelet transform (scaling contribution in real space).
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
 * Perform wavelet transform in real space (from scratch, gives pixel space components).
 * Input function is real.
 * Sampling scheme : MW sampling.
 * Spherical wavelets : analysis in real space.
 *
 * \param[out]  f_wav Wavelet transform (wavelet contribution in real space).
 * \param[out]  f_scal Wavelet transform (scaling contribution in real space).
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
 * Perform wavelet transform in real space (from scratch, gives pixel space components).
 * Input function is real.
 * Sampling scheme : MW sampling.
 * Spherical wavelets : synthesis in real space.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Wavelet transform (wavelet contribution in real space).
 * \param[in]  f_scal Wavelet transform (scaling contribution in real space).
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

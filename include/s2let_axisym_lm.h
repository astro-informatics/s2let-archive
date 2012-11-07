// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_AXISYM_LM
#define S2LET_AXISYM_LM

#include <complex.h> 


/*!
 * Allocates arrays for wavelets and scaling functions in harmonic space.
 *
 * \param[out]  f_wav_lm Harmonic coefficients of the wavelets. Each wavelet has size L*L and there are (J+1-J_min) scales.
 * \param[out]  f_scal_lm Harmonic coefficients of the scaling function (L*L).
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_lm_allocate_f_wav(complex double **f_wav_lm, complex double **f_scal_lm, int B, int L, int J_min);

/*!
 * Allocates multiresolution arrays for wavelets and scaling functions in harmonic space.
 *
 * \param[out]  f_wav_lm Harmonic coefficients of the wavelets. The size of each wavelet depends on its band-limit.
 * \param[out]  f_scal_lm Harmonic coefficients of the scaling function (L*L).
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_lm_allocate_f_wav_multires(complex double **f_wav_lm, complex double **f_scal_lm, int B, int L, int J_min);


/*!
 * Allocates arrays for the kernels of the wavelets and the scaling functions.
 *
 * \param[out]  wav_lm Wavelet kernels.
 * \param[out]  scal_lm Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_axisym_lm_allocate_wav(double **wav_lm, double **scal_lm, int B, int L);

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
 void s2let_axisym_lm_wav(double *wav_lm, double *scal_lm, int B, int L, int J_min);


/*!
 * Spherical wavelets : full resolution analysis in harmonic space.
 * Perform wavelet transform from precomputed kernels and gives the harmonic coefficients.
 *
 * \param[out]  f_wav_lm Wavelet transform (harmonic coefficients of wavelet contribution).
 * \param[out]  f_scal_lm Wavelet transform (harmonic coefficients of scaling contribution).
 * \param[in]  flm Spherical harmonic decomposition of input function.
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_lm Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_lm_wav_analysis(complex double *f_wav_lm, complex double *f_scal_lm, const complex double *flm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min);

/*!
 * Spherical wavelets : full resolution synthesis in harmonic space.
 * Perform wavelet transform in harmonic space from precomputed kernels and gives harmonic coefficients.
 *
 * \param[out]  flm Spherical harmonic decomposition of input function.
 * \param[in]  f_wav_lm Wavelet transform (harmonic coefficients of wavelet contribution).
 * \param[in]  f_scal_lm Wavelet transform (harmonic coefficients of scaling contribution).
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_lm Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
 void s2let_axisym_lm_wav_synthesis(complex double *flm, const complex double *f_wav_lm, const complex double *f_scal_lm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min);

/*!
 * Spherical wavelets : multiresolution analysis in harmonic space.
 * Perform multiresolution wavelet transform in harmonic space from precomputed kernels and gives harmonic coefficients.
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
void s2let_axisym_lm_wav_analysis_multires(complex double *f_wav_lm, complex double *f_scal_lm, const complex double *flm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min);

/*!
 * Spherical wavelets : multiresolution synthesis in harmonic space.
 * Perform multiresolution wavelet transform in harmonic space from precomputed kernels and gives harmonic coefficients.
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
 void s2let_axisym_lm_wav_synthesis_multires(complex double *flm, const complex double *f_wav_lm, const complex double *f_scal_lm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min);

#endif

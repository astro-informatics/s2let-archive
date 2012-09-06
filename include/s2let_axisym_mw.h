// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_AXISYM_MW
#define S2LET_AXISYM_MW

#include <complex.h> 

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
void s2let_axisym_allocate_f_wav_multires(complex double **f_wav, complex double **f_scal, int B, int L, int J_min);

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
void s2let_axisym_allocate_f_wav_multires_real(double **f_wav, double **f_scal, int B, int L, int J_min);

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
void s2let_axisym_allocate_f_wav(complex double **f_wav, complex double **f_scal, int B, int L, int J_min);

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
void s2let_axisym_allocate_f_wav_real(double **f_wav, double **f_scal, int B, int L, int J_min);

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
void s2let_axisym_wav_analysis(complex double *f_wav, complex double *f_scal, const complex double *f, int B, int L, int J_min);

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
void s2let_axisym_wav_synthesis(complex double *f, const complex double *f_wav, const complex double *f_scal, int B, int L, int J_min);

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
void s2let_axisym_wav_analysis_real(double *f_wav, double *f_scal, const double *f, int B, int L, int J_min);

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
void s2let_axisym_wav_synthesis_real(double *f, const double *f_wav, const double *f_scal, int B, int L, int J_min);

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
void s2let_axisym_wav_analysis_multires(complex double *f_wav, complex double *f_scal, const complex double *f, int B, int L, int J_min);

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
void s2let_axisym_wav_synthesis_multires(complex double *f, const complex double *f_wav, const complex double *f_scal, int B, int L, int J_min);

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
void s2let_axisym_wav_analysis_multires_real(double *f_wav, double *f_scal, const double *f, int B, int L, int J_min);

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
void s2let_axisym_wav_synthesis_multires_real(double *f, const double *f_wav, const double *f_scal, int B, int L, int J_min);

/*!
 * Treshold real wavelets in real space, MW sampling, multiresolution.
 *
 * \param[inout]  g_wav Array of wavelets maps, MW sampling.
 * \param[in]  treshold A treshold rule, i.e. a number for every scale j.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_wav_hardthreshold_multires_real(double *g_wav, const double *treshold, int B, int L, int J_min);

/*!
 * Treshold real wavelets in real space, MW sampling, full resolution.
 *
 * \param[inout]  g_wav Array of wavelets maps, MW sampling.
 * \param[in]  treshold A treshold rule, i.e. a number for every scale j.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_axisym_wav_hardthreshold_real(double *g_wav, const double *treshold, int B, int L, int J_min);

#endif

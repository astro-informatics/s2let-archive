// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_TILING
#define S2LET_TILING

#define PI    3.141592653589793238462643383279502884197

/*!
 * Computes band-limit of a specific wavelet scale.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  j Wavelet scale.
 * \retval band-limit
 */
int s2let_bandlimit(int B, int j);

/*!
 * Computes minimum harmonic index supported by needlets.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval ell_min
 */
int s2let_el_min(int B, int J_min);

/*!
 * Computes needlet maximum level required to ensure exact reconstruction.
 *
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  B Wavelet parameter.
 * \retval j_max
 */
int s2let_j_max(int L, int B);

/*!
 * Allocates tiling in harmonic space.
 *
 * \param[out]  kappa Kernel functions for the wavelets.
 * \param[out]  kappa0 Kernel for the scaling function.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_tiling_axisym_allocate(double **kappa, double **kappa0, int B, int L);

/*!
 * Generates tiling in harmonic space.
 *
 * \param[out]  kappa Kernel functions for the wavelets.
 * \param[out]  kappa0 Kernel for the scaling function.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_tiling_axisym(double *kappa, double *kappa0, int B, int L, int J_min);

/*!
 * Generates smooth functions to construct the tiling.
 *
 * \param[out]  phi2 Smooth step functions for the wavelets.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_tiling_phi2(double *phi2, int B, int L, int J_min);

/*!
 * Checks exactness of the harmonic tiling.
 *
 * \param[in]  kappa Kernel functions for the wavelets.
 * \param[in]  kappa0 Kernel for the scaling function.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval Achieved accuracy (should be lower than e-12).
 */
double s2let_tiling_axisym_check_identity(double *kappa, double *kappa0, int B, int L, int J_min);

#endif

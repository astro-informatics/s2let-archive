// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_TILING
#define S2LET_TILING

#define PI    3.141592653589793238462643383279502884197

void s2let_switch_wavtype(int typenum);

int s2let_bandlimit(int j, const s2let_parameters_t *parameters);

int s2let_el_min(int B, int J_min);

int s2let_j_max(int L, int B);

void s2let_tiling_axisym_allocate(double **kappa, double **kappa0, int B, int L);

void s2let_tiling_axisym(double *kappa, double *kappa0, int B, int L, int J_min);

void s2let_tiling_direction_allocate(complex double **s_elm, int L, int N);

void s2let_tiling_direction(complex double *s_elm, int L, int N);

void s2let_tiling_wavelet_allocate(complex double **psi, double **phi, int B, int L, int N);

void s2let_tiling_wavelet(
    complex double *psi,
    double *phi,
    int B,
    int L,
    int J_min,
    int N,
    int spin,
    s2let_wav_norm_t normalisation,
    int original_spin
);

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

double s2let_tiling_axisym_check_identity(double *kappa, double *kappa0, int B, int L, int J_min);

double s2let_tiling_direction_check_identity(complex double *s_elm, int L, int N);

double s2let_tiling_wavelet_check_identity(complex double *psi, double *phi, int B, int L, int J_min, int N, int spin);

#endif

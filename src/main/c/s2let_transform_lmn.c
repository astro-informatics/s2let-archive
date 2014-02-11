// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

static inline int lm2ind(int el, int m)
{
    return el*el + el + m;
}

static inline int lmn2ind(int el, int m, int n, int L, int N)
{
    return (N-1+n) * L*L + el*el + el + m;
}

/*!
 * Allocates arrays for directional wavelet transform in harmonic space.
 *
 * \param[out]  f_wav_lmn Wigner coefficients of the wavelet contributions.
 *                        Each wavelet has size (2*N-1)*L*L and there are
 *                        (J-J_min+1) scales.
 * \param[out]  f_scal_lm Spherical harmonic coefficients of the scaling
 *                        contribution (L*L).
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \retval none
 */
void s2let_allocate_f_wav_lmn(
    complex double **f_wav_lmn,
    complex double **f_scal_lm,
    int B,
    int L,
    int J_min,
    int N
) {
    int J = s2let_j_max(L, B);
    *f_wav_lmn = calloc((J-J_min+1) * (2*N-1) * L*L, sizeof **f_wav_lmn);
    *f_scal_lm = calloc(L*L, sizeof **f_scal_lm);
}

/*!
 * Allocates multi-resolution arrays for directional wavelet transform in
 * harmonic space.
 *
 * \param[out]  f_wav_lmn Wigner coefficients of the wavelet contributions.
 *                        The size of each wavelet depends on its band-limit.
 * \param[out]  f_scal_lm Spherical harmonic coefficients of the scaling
 *                        contribution (L*L).
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \retval none
 */
void s2let_allocate_f_wav_lmn_multires(
    complex double **f_wav_lmn,
    complex double **f_scal_lm,
    int B,
    int L,
    int J_min,
    int N
) {
    int J = s2let_j_max(L, B);
    int j, bandlimit, total = 0;
    for (j = J_min; j <= J; ++j)
    {
        bandlimit = MIN(s2let_bandlimit(j, J_min, B, L), L);
        total += (2*N-1) * bandlimit * bandlimit;
    }
    *f_wav_lmn = calloc(total, sizeof **f_wav_lmn);
    bandlimit = MIN(s2let_bandlimit(J_min-1, J_min, B, L), L);
    *f_scal_lm = calloc(bandlimit * bandlimit, sizeof **f_scal_lm);
}

/*!
 * Spherical wavelets: full resolution analysis in harmonic space.
 * Perform directional wavelet transform from precomputed kernels
 * to give Wigner coefficients.
 *
 * \param[out]  f_wav_lmn Wavelet transform (Wigner coefficients of wavelet contribution).
 * \param[out]  f_scal_lm Wavelet transform (spherical harmonic coefficients of scaling contribution).
 * \param[in]  flm Spherical harmonic decomposition of input function.
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_l Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_wav_analysis_harmonic(
    complex double *f_wav_lmn,
    complex double *f_scal_lm,
    const complex double *flm,
    const complex double *wav_lm,
    const double *scal_l,
    int B,
    int L,
    int J_min,
    int N
) {
    int j, el, m ,n;
    int J = s2let_j_max(L, B);

    complex double psi;
    double phi;

    int offset = 0;

    for (j = J_min; j <= J; j++)
    {
        for (n = -N+1; n < N; ++n)
        {
            for (el = ABS(n); el < L; ++el)
            {
                psi = conj(wav_lm[j*L*L + el*el + el + n]);
                for (m = -el; m <= el; ++m)
                {
                    f_wav_lmn[offset + lmn2ind(el,m,n,L,N)] =
                        8*PI*PI/(2*el+1) * flm[lm2ind(el,m)] * psi;
                }
            }
        }
        offset += (2*N-1) * L*L;
    }

    for (el = 0; el < L; ++el)
    {
        phi = scal_l[el];
        for (m = -el; m <= el; ++m)
        {
            f_scal_lm[lm2ind(el,m)] =
                flm[lm2ind(el,m)] * phi;
        }
    }
}


/*!
 * Spherical wavelets: full resolution synthesis in harmonic space.
 * Perform directional wavelet transform from precomputed kernels
 * to give Wigner coefficients.
 *
 * \param[out]  flm Spherical harmonic decomposition of input function.
 * \param[in]  f_wav_lmn Wavelet transform (Wigner coefficients of wavelet contribution).
 * \param[in]  f_scal_lm Wavelet transform (spherical harmonic coefficients of scaling contribution).
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_l Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_wav_synthesis_harmonic(
    complex double *flm,
    const complex double *f_wav_lmn,
    const complex double *f_scal_lm,
    const complex double *wav_lm,
    const double *scal_l,
    int B,
    int L,
    int J_min,
    int N
) {
    int j, el, m ,n;
    int J = s2let_j_max(L, B);

    complex double psi;
    double phi;

    int offset = 0;


    for (j = J_min; j <= J; ++j)
    {
        for (n = -N+1; n < N; ++n)
        {
            for (el = ABS(n); el < L; ++el)
            {
                psi = wav_lm[j*L*L + el*el + el + n];
                for (m = -el; m <= el; ++m)
                {
                    flm[lm2ind(el,m)] +=
                        (2*el+1)/(8*PI*PI) *
                        f_wav_lmn[offset + lmn2ind(el,m,n,L,N)] * psi;
                }
            }
        }
        offset += (2*N-1) * L*L;
    }

    for (el = 0; el < L; ++el)
    {
        phi = scal_l[el];
        for (m = -el; m <= el; ++m)
        {
            flm[lm2ind(el,m)] +=
                f_scal_lm[lm2ind(el,m)] * phi;
        }
    }
}

/*!
 * Spherical wavelets: multi-resolution analysis in harmonic space.
 * Perform directional wavelet transform from precomputed kernels
 * to give Wigner coefficients.
 *
 * \param[out]  f_wav_lmn Wavelet transform (Wigner coefficients of wavelet contribution).
 * \param[out]  f_scal_lm Wavelet transform (spherical harmonic coefficients of scaling contribution).
 * \param[in]  flm Spherical harmonic decomposition of input function.
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_l Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_wav_analysis_harmonic_multires(
    complex double *f_wav_lmn,
    complex double *f_scal_lm,
    const complex double *flm,
    const complex double *wav_lm,
    const double *scal_l,
    int B,
    int L,
    int J_min,
    int N
) {
    int j, el, m ,n;
    int J = s2let_j_max(L, B);
    int bandlimit;

    complex double psi;
    double phi;

    int offset = 0;

    for (j = J_min; j <= J; j++)
    {
        bandlimit = MIN(s2let_bandlimit(j, J_min, B, L), L);
        for (n = -N+1; n < N; ++n)
        {
            for (el = ABS(n); el < bandlimit; ++el)
            {
                psi = conj(wav_lm[j*L*L + el*el + el + n]);
                for (m = -el; m <= el; ++m)
                {
                    f_wav_lmn[offset + lmn2ind(el,m,n,bandlimit,N)] =
                        8*PI*PI/(2*el+1) * flm[lm2ind(el,m)] * psi;
                }
            }
        }
        offset += (2*N-1) * bandlimit*bandlimit;
    }

    bandlimit = MIN(s2let_bandlimit(J_min-1, J_min, B, L), L);
    for (el = 0; el < bandlimit; ++el)
    {
        phi = scal_l[el];
        for (m = -el; m <= el; ++m)
        {
            f_scal_lm[lm2ind(el,m)] =
                flm[lm2ind(el,m)] * phi;
        }
    }
}


/*!
 * Spherical wavelets: multi-resolution synthesis in harmonic space.
 * Perform directional wavelet transform from precomputed kernels
 * to give Wigner coefficients.
 *
 * \param[out]  flm Spherical harmonic decomposition of input function.
 * \param[in]  f_wav_lmn Wavelet transform (Wigner coefficients of wavelet contribution).
 * \param[in]  f_scal_lm Wavelet transform (spherical harmonic coefficients of scaling contribution).
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_l Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_wav_synthesis_harmonic_multires(
    complex double *flm,
    const complex double *f_wav_lmn,
    const complex double *f_scal_lm,
    const complex double *wav_lm,
    const double *scal_l,
    int B,
    int L,
    int J_min,
    int N
) {
    int j, el, m ,n;
    int J = s2let_j_max(L, B);
    int bandlimit;

    complex double psi;
    double phi;

    int offset = 0;


    for (j = J_min; j <= J; ++j)
    {
        bandlimit = MIN(s2let_bandlimit(j, J_min, B, L), L);
        for (n = -N+1; n < N; ++n)
        {
            for (el = ABS(n); el < bandlimit; ++el)
            {
                psi = wav_lm[j*L*L + el*el + el + n];
                for (m = -el; m <= el; ++m)
                {
                    flm[lm2ind(el,m)] +=
                        (2*el+1)/(8*PI*PI) *
                        f_wav_lmn[offset + lmn2ind(el,m,n,bandlimit,N)] * psi;
                }
            }
        }
        offset += (2*N-1) * bandlimit*bandlimit;
    }

    bandlimit = MIN(s2let_bandlimit(J_min-1, J_min, B, L), L);
    for (el = 0; el < bandlimit; ++el)
    {
        phi = scal_l[el];
        for (m = -el; m <= el; ++m)
        {
            flm[lm2ind(el,m)] +=
                f_scal_lm[lm2ind(el,m)] * phi;
        }
    }
}

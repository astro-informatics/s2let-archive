// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_HPX
#define S2LET_HPX

#include <complex.h>

/*!
 * Restore real healpix map from spherical harmonic coefficients.
 * Interface for HEALPIX Fortran code alm2map.
 *
 * \param[out]  f Output healpix map.
 * \param[in]  flm Spherical harmonic coefficients.
 * \param[in]  nside Healpix resolution of the output map.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_hpx_alm2map_real(double* f, const complex double* flm, int nside, int L);


/*!
 * Compute spherical harmonic coefficients of a real healpix map.
 * Interface for HEALPIX Fortran code map2alm.
 *
 * \param[out]  flm Spherical harmonic coefficients.
 * \param[in]  f Input healpix map.
 * \param[in]  nside Healpix resolution of the output map.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_hpx_map2alm_real(complex double* flm, const double* f, int nside, int L);


/*!
 * Read Healpix map from a FITS file.
 * Interface for HEALPIX Fortran code read_bintab.
 *
 * \param[out]  f Input healpix map.
 * \param[in]  file Filename.
 * \param[in]  nside Healpix resolution of the output map.
 * \retval none
 */
void s2let_hpx_read_map(double* f, char* file, int nside);


/*!
 * Write Healpix map to a FITS file.
 * Interface for HEALPIX Fortran code write_bintab.
 *
 * \param[in]  file Filename.
 * \param[in]  f Input healpix map.
 * \param[in]  nside Healpix resolution of the output map.
 * \retval none
 */
void s2let_hpx_write_map(char* file, const double* f, int nside);


/*!
 * Allocate Healpix map.
 *
 * \param[inout]  f healpix map.
 * \param[in]  nside Healpix resolution.
 * \retval none
 */
void s2let_hpx_allocate_real(double **f, int nside);

#endif

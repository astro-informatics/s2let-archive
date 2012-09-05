// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_HPXTOOLS
#define S2LET_HPXTOOLS

#include <complex.h>

// Fortran interfaces to Healpix F90 library ; see s2let_hpx.f90
extern void healpix_inverse_real_();
extern void healpix_forward_real_();
extern void s2let_hpx_write_map_();
extern void s2let_hpx_read_map_();

void s2let_hpx_alm2map_real(double* f, const complex double* flm, int nside, int L);
void s2let_hpx_map2alm_real(complex double* flm, const double* f, int nside, int L);

void read_healpix_map(double* f, char* file, int nside);
void write_healpix_map(char* file, const double* f, int nside);

#endif

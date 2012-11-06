// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <stdlib.h>
#include <complex.h> 

// Fortran interfaces to Healpix F90 library ; see s2let_hpx.f90
extern void healpix_inverse_real_();
extern void healpix_forward_real_();
extern void write_healpix_map_();
extern void read_healpix_map_();

void s2let_hpx_alm2map_real(double* f, const complex double* flm, int nside, int L)
{
  healpix_inverse_real_(f, flm, &nside, &L);
}

void s2let_hpx_map2alm_real(complex double* flm, const double* f, int nside, int L)
{
  healpix_forward_real_(flm, f, &nside, &L);
}

void s2let_read_hpx_map(double* f, char* file, int nside)
{
  s2let_check_hpx_ordering(file);
  read_healpix_map_(f, file, &nside);
}

void s2let_write_hpx_map(char* file, const double* f, int nside)
{
  write_healpix_map_(file, f, &nside);
}

void s2let_allocate_hpx_real(double **f, int nside)
{
  *f = (double*)calloc(12*nside*nside, sizeof(double));
}

// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <stdlib.h>
#include <complex.h>

// Fortran interfaces to Healpix F90 library ; see s2let_hpx.f90
extern void healpix_inverse_real_(double *, const complex double *, int *, int *);
extern void healpix_forward_real_(complex double *, const double *, int *, int *);
extern void healpix_inverse_spin_real_(double *, double *, const complex double *, const complex double *, int *, int *, int *);
extern void healpix_forward_spin_real_(complex double *, complex double *, const double *,const double *,int *, int *, int *);
extern void write_healpix_map_(char *, const double *, int *);
extern void read_healpix_map_(double *, char *, int *);
extern void read_healpix_maps_(double *, char *, int *, int*);

void s2let_hpx_alm2map_real(double* f, const complex double* flm, int nside, int L)
{
  healpix_inverse_real_(f, flm, &nside, &L);
}

void s2let_hpx_map2alm_real(complex double* flm, const double* f, int nside, int L)
{
  healpix_forward_real_(flm, f, &nside, &L);
}

void s2let_hpx_alm2map_spin_real(double* fQ, double* fU, const complex double* flmE, const complex double* flmB, int nside, int L, int spin)
{
  healpix_inverse_spin_real_(fQ, fU, flmE, flmB, &nside, &L, &spin);
}

void s2let_hpx_map2alm_spin_real(complex double* flmE, complex double* flmB, const double* fQ, const double* fU, int nside, int L, int spin)
{
  healpix_forward_spin_real_(flmE, flmB, fQ, fU, &nside, &L, &spin);
}

void s2let_hpx_read_maps(double* f, char* file, int nside, int nmaps)
{
  read_healpix_maps_(f, file, &nside, &nmaps);
}

void s2let_hpx_read_map(double* f, char* file, int nside)
{
  read_healpix_map_(f, file, &nside);
}

void s2let_hpx_write_map(char* file, const double* f, int nside)
{
  write_healpix_map_(file, f, &nside);
}

void s2let_hpx_allocate_real(double **f, int nside)
{
  *f = calloc(12*nside*nside, sizeof **f);
}

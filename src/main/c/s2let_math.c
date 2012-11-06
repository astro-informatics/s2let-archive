// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*!
 * Tiling function for S2DW wavelets
 */
double f_s2dw(double k, int B)
{
  double t = (k - (1 / (double)B)) * (2.0 * B / (double)(B-1)) - 1;
  return exp(-2.0 / (1.0 - pow(t, 2.0))) / k;
}

/*!
 * Tiling function for needlets
 */
double f_needlet(double t)
{
  return exp(-1.0 / (1.0 - pow(t, 2.0))) ;
}


/*!
 * Computes smooth "Schwarts" functions.
 */
double s2let_kappa0_quadtrap_s2dw(double a, double b, int n, int B)
{
  double sum = 0;
  double f1, f2;
  int i;
  double h = (b - a) / n;

  if( a == b ){
    return 0;
  }else{
    for (i = 0; i < n; i++){
      f1 = f_s2dw(a + i * h, B);
      f2 = f_s2dw(a + (i + 1) * h, B);
      if(!isnan(f1) && !isinf(f1) && !isnan(f2) && !isinf(f2))
	sum += ((f1 + f2) * h) / 2;
    }
  }
  return sum;
}

/*!
 * Computes smooth "Schwarts" functions.
 */
double s2let_kappa0_quadtrap_needlet(double a, double b, int n)
{
  double sum = 0;
  double f1, f2;
  int i;
  double h = (b - a) / n;

  if( a == b ){
    return 0;
  }else{
    for (i = 0; i < n; i++){
      f1 = f_needlet(a + i * h);
      f2 = f_needlet(a + (i + 1) * h);
      if(!isnan(f1) && !isinf(f1) && !isnan(f2) && !isinf(f2))
	sum += ((f1 + f2) * h) / 2;
    }
  }
  return sum;
}


/*!
 * Random number from seed (Numerical Recipes)
 */
double ran2_dp(int idum) {
  int IM1=2147483563,IM2=2147483399,IMM1=IM1-1, 
    IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, 
    NTAB=32,NDIV=1+IMM1/NTAB;

  double AM=1./IM1,EPS=1.2e-7,RNMX=1.-EPS;
  int j,k;
  static int iv[32],iy,idum2 = 123456789; 
  // N.B. in C static variables are initialised to 0 by default.

  if (idum <= 0) {
    idum= (-idum>1 ? -idum : 1); // max(-idum,1);
    idum2=idum;
    for(j=NTAB+8;j>=1;j--) {
      k=idum/IQ1;
      idum=IA1*(idum-k*IQ1)-k*IR1;
      if (idum < 0) idum=idum+IM1;
      if (j < NTAB) iv[j-1]=idum;
    }
    iy=iv[0];
  }
  k=idum/IQ1;
  idum=IA1*(idum-k*IQ1)-k*IR1;
  if (idum < 0) idum=idum+IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2=idum2+IM2;
  j=1+iy/NDIV;
  iy=iv[j-1]-idum2;
  iv[j-1]=idum;
  if(iy < 1)iy=iy+IMM1;
  return (AM*iy < RNMX ? AM*iy : RNMX); // min(AM*iy,RNMX);
}


/*!
 * Generate random harmonic coefficients for a complex map.
 *
 * \param[out]  flm Harmonic coefficients.
 * \param[in]  L Band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_axisym_random_flm(complex double *flm, int L, int seed)
{
  int i;
  srand( time(NULL) );
  for (i=0; i<L*L; i++){
    flm[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
  }
}

/*!
 * Generate random harmonic coefficients corresponding to a real map.
 *
 * \param[out]  flm Harmonic coefficients.
 * \param[in]  L Band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_axisym_random_flm_real(complex double *flm, int L, int seed) {
  int el, m, msign, i, i_op;
  for (el=0; el<L; el++) {
    m = 0;
    i = el*el + el + m ;
    flm[i] = (2.0*ran2_dp(seed) - 1.0);
    for (m=1; m<=el; m++) {
      i = el*el + el + m ;
      flm[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
      i_op = el*el + el - m ;
      msign = m & 1;
      msign = 1 - msign - msign; // (-1)^m
      flm[i_op] = msign * conj(flm[i]);
    }
  }
}

/*!
 * Max absolute error between two complex arrays
 */
double maxerr_cplx(complex double *a, complex double *b, int size)
{
  double value = 0;
  int i;
  for(i = 0; i<size; i++){
    value = MAX( cabs( a[i]-b[i] ), value );
  }
  return value;
}

/*!
 * Max absolute error between two real arrays
 */
double maxerr(double *a, double *b, int size)
{
  double value = 0;
  int i;
  for(i = 0; i<size; i++){
    value = MAX( abs( a[i]-b[i] ), value );
  }
  return value;
}

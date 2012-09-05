// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <math.h>

/*!
 * Generating function for the smooth "Schwarts" functions.
 *
 */
double f(double k, int B)
{
	double t = (k - (1 / (double)B)) * (2.0 * B / (double)(B-1)) - 1;
    return exp(-2.0 / (1.0 - pow(t, 2.0))) / k;
}


/*!
 * Computes smooth "Schwarts" functions.
 *
 */
double s2let_kappa0_quadtrap(double a, double b, int n, int B)
{
    double sum = 0;
    double f1, f2;
    int i;
    double h = (b - a) / n;

    if( a == b ){
    	return 0;
    }else{
	    for (i = 0; i < n; i++){
	    	f1 = f(a + i * h, B);
	    	f2 = f(a + (i + 1) * h, B);
        if(!isnan(f1) && !isinf(f1) && !isnan(f2) && !isinf(f2))
          sum += ((f1 + f2) * h) / 2;
	    }
	}
  return sum;
}

/*!
 * Simpson rule for numerical integration.
 *
 */
double simpson(double a, double b, int n, int B)
{
  long double integral,x,h;
  long double part,coeff;
  int i;
  part = (b-a) / (long double) n;
  h = part / (long double) 3.0;

  integral = 0;
  x = 0;

  for (i=0;i<n;i++)
  {
     integral = integral+f(x,B)+3.0*f(x+h,B)+3.0*f(x+2*h,B)+f(x+3.0*h,B);
     x = x + 3.0*h ;
  }
  coeff = 3.0*h / 8.0 ;
  integral = coeff * integral ;
  return integral;
}

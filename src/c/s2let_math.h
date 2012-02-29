// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_MATH
#define S2LET_MATH

#define PI    3.141592653589793238462643383279502884197

double s2let_kappa0_integrand(double k, int B);

double s2let_kappa0_quadtrap(double a, double b, int n, int B);

#endif
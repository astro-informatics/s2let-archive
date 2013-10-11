// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <math.h>
#include <stdlib.h>

//
//
typedef enum { S2DW, NEEDLET, SPLINE } s2let_kernel_type ;
s2let_kernel_type s2let_kernel = S2DW;
//
//

void s2let_switch_wavtype(int typenum)
{
  //printf("Input wavelet type : %i\n",typenum);
  if (typenum == 1){
    //printf("Kernel switch 1: using scale-discretised wavelets.\n");
    s2let_kernel = S2DW;
  } else if (typenum == 2){
    //printf("Kernel switch 2: using needlets.\n");
    s2let_kernel = NEEDLET;
  } else if (typenum == 3){
    //printf("Kernel switch 3: using cubic splines wavelets.\n");
    s2let_kernel = SPLINE;
  } else {
    printf("Kernel number should be 1, 2 or 3. Default is 1.\n");
    s2let_kernel = S2DW;
  }
}

int s2let_bandlimit(int j, int J_min, int B, int L)
{
  if(s2let_kernel == S2DW || s2let_kernel == NEEDLET)
    return ceil(pow(B, j+1));
  if(s2let_kernel == SPLINE){
    int Jmax = s2let_j_max(L, B);
    if (j == Jmax) return L;
    //if (j < J_min) return ceil(L / (double) pow(B, Jmax-J_min-1));
    return ceil(L / (double) pow(B, Jmax-j-2)); 
  }
}

int s2let_el_min(int B, int J_min)
{
  if(s2let_kernel == S2DW || s2let_kernel == NEEDLET)
    return ceil(pow(B, J_min));
  if(s2let_kernel == SPLINE)
    return 0;
}

int s2let_j_max(int L, int B)
{
  return ceil(log(L) / log(B));
}

void s2let_tiling_axisym_allocate(double **kappa, double **kappa0, int B, int L)
{
  int J = s2let_j_max(L, B);
  *kappa = (double*)calloc((J+1) * L, sizeof(double));
  *kappa0 = (double*)calloc(L, sizeof(double));
}

void s2let_tiling_phi2_s2dw(double *phi2, int B, int L, int J_min)
{
  int j, l;
  int J = s2let_j_max(L, B);
  int n = 300;

  double kappanorm = s2let_math_kappa0_quadtrap_s2dw(1.0 / (double)B, 1.0, n, B);
  for (j = 0; j <= J+1; j++){
    for (l = 0; l < L; l++){
      if (l < pow(B,j-1)) {
        phi2[l+j*L] = 1;
      } else if (l > pow(B,j)) {
        phi2[l+j*L] = 0;
      } else {
	      phi2[l+j*L] = s2let_math_kappa0_quadtrap_s2dw((double)l / pow(B, j), 1.0, n, B) / kappanorm;
      }
    }
  }
}

void s2let_tiling_phi2_needlet(double *phi2, int B, int L, int J_min)
{
  int j, l;
  int J = s2let_j_max(L, B);
  int n = 300;
  double u;

  double kappanorm = s2let_math_kappa0_quadtrap_needlet(-1.0, 1.0, n);
  for (j = 0; j <= J+1; j++){
    for (l = 0; l < L; l++){
      if (l < pow(B,j-1)) {
	       phi2[l+j*L] = 1;
      } else if (l > pow(B,j)) {
	       phi2[l+j*L] = 0;
      } else {
	       u = 1.0 - 2.0 * B / (B - 1.0) * ( l * pow(B, -j) - 1.0 / B );
	       phi2[l+j*L] = s2let_math_kappa0_quadtrap_needlet(-1.0, u, n) / kappanorm;
      }
    }
  }
}

void s2let_tiling_phi2_spline(double *phi2, int B, int L, int J_min)
{
  int j = 0, l;
  int J = s2let_j_max(L, B);
  phi2[(J+1-j)*L] = 1.0;
  for (l = 1; l < L; l++){
      phi2[l+(J+1-j)*L] = 1.0;
    }
  for (j = 1; j <= J+1; j++){
    double bl = (double) L / (double) pow(B, j-2);
    phi2[(J+1-j)*L] = 1.0;
    for (l = 1; l < L; l++){
        if (l > bl)
          phi2[l+(J+1-j)*L] = 0.0;
        else
          phi2[l+(J+1-j)*L] = s2let_math_spline_scalingfct((double) l, bl);
      }
  }
}

void s2let_tiling_axisym(double *kappa, double *kappa0, int B, int L, int J_min)
{
  int j, l;
  int J = s2let_j_max(L, B);

  double previoustemp = 0.0, temp;
  double *phi2 = (double*)calloc((J+2) * L, sizeof(double));

  if(s2let_kernel == SPLINE)
    s2let_tiling_phi2_spline(phi2, B, L, J_min); // SPLINE tiling
  if(s2let_kernel == S2DW)
    s2let_tiling_phi2_s2dw(phi2, B, L, J_min); // S2DW tiling
  if(s2let_kernel == NEEDLET)
    s2let_tiling_phi2_needlet(phi2, B, L, J_min); // Needlet tiling

  for (l = 0; l < L; l++)
      kappa0[l] = sqrt(phi2[l+J_min*L]);
  
  for (j = J_min; j <= J; j++){
    for (l = 0; l < L; l++){
      temp = sqrt(phi2[l+(j+1)*L] - phi2[l+j*L]);
      if( isnan(temp) || isinf(temp) )
        kappa[l+j*L] = previoustemp;
      else
        kappa[l+j*L] = temp;
      previoustemp = temp;
    }
    for (l = 0; l < L; l++)
      if( !finite(kappa[l+j*L]) ) 
        kappa[l+j*L] = kappa[l+j*L-1];
  }
  free(phi2);
}

double s2let_tiling_axisym_check_identity(double *kappa, double *kappa0, int B, int L, int J_min)
{
  int l, j;
  int J = s2let_j_max(L, B);
  //int l_min = s2let_el_min(B, J_min);
  double sum = 0;
		
  double *ident;
  ident = (double*)calloc(L, sizeof(double));

  //printf("Tilling identity: ");
  for (l = 0; l < L; l++){
    ident[l] = pow(kappa0[l], 2.0);
    //sum += ident[l] - 1.0;
  }
  //printf(" %2.2f ", ident[0]);
  for (l = 0; l < L; l++){
    for (j = J_min; j <= J; j++){
      ident[l] += pow(kappa[l+j*L], 2.0);
    }
    //printf(" %2.2f ", ident[l]);
    sum += ident[l] - 1.0;
  }

  return sum;
  free(ident);
}


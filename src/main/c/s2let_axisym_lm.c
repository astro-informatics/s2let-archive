// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <complex.h> 
#include <math.h>
#include <stdlib.h>

int lm2ind(int el, int em)
{
  return el*el + el + em;
}

void s2let_axisym_lm_allocate_f_wav(complex double **f_wav_lm, complex double **f_scal_lm, int B, int L, int J_min)
{
  int J = s2let_j_max(L, B);
  *f_wav_lm = (complex double*)calloc((J+1-J_min) * L * L, sizeof(complex double));
  *f_scal_lm = (complex double*)calloc(L * L, sizeof(complex double));
}

void s2let_axisym_lm_allocate_f_wav_multires(complex double **f_wav_lm, complex double **f_scal_lm, int B, int L, int J_min)
{
  int J = s2let_j_max(L, B);
  int j, bandlimit, total = 0;
  for(j = J_min; j <= J; j++){
    bandlimit = MIN(s2let_bandlimit(j, J_min, B, L), L);
    total += bandlimit * bandlimit;
  }
  *f_wav_lm = (complex double*)calloc(total, sizeof(complex double));
  bandlimit = MIN(s2let_bandlimit(J_min-1, J_min, B, L), L);
  *f_scal_lm = (complex double*)calloc(bandlimit * bandlimit, sizeof(complex double));
}

void s2let_axisym_lm_allocate_wav(double **wav_lm, double **scal_lm, int B, int L)
{
  int J = s2let_j_max(L, B);
  *wav_lm = (double*)calloc((J+1) * L, sizeof(double));
  *scal_lm = (double*)calloc(L, sizeof(double));
}

void s2let_axisym_lm_wav(double *wav_lm, double *scal_lm, int B, int L, int J_min)
{
  int j, l;
  int J = s2let_j_max(L, B);
  //int J_min = 0;
  //int l_min = s2let_axisym_el_min(B, J_min);
  double k0;
  double *kappa, *kappa0;
  s2let_tiling_axisym_allocate(&kappa, &kappa0, B, L);
  s2let_tiling_axisym(kappa, kappa0, B, L, J_min);

  for (j = J_min; j <= J; j++){
    for (l = 0; l < L; l++){
      k0 = sqrt( (2 * l + 1) / (4.0 * PI) ) * kappa[l+j*L];
      wav_lm[j*L+l] = k0;
    }
  }
  for (l = 0; l < L; l++){
    k0 = sqrt( (2 * l + 1) / (4.0 * PI) ) * kappa0[l];
    scal_lm[l] = k0;
  }

  free(kappa);
  free(kappa0);
}

void s2let_axisym_lm_wav_analysis(complex double *f_wav_lm, complex double *f_scal_lm, const complex double *flm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min)
{
  int offset, j, l, m;
  int J = s2let_j_max(L, B);
  double wav0, scal0;
  //int l_min = s2let_axisym_el_min(B, J_min);

  offset = 0;
  for (j = J_min; j <= J; j++){
    for (l = 0; l < L; l++){
      wav0 = sqrt((4.0*PI)/(2.0*l+1.0)) * wav_lm[j*L+l];
      for (m = -l; m <= l; m++){
	      f_wav_lm[offset + l*l + l + m] = flm[lm2ind(l,m)] * wav0 ;
      }
    }
    offset += L * L;
  }
  for (l = 0; l < L; l++){
    scal0 = sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lm[l];
    for (m = -l; m <= l; m++){
      f_scal_lm[lm2ind(l,m)] = flm[lm2ind(l,m)] * scal0 ;
    }
  }
}

void s2let_axisym_lm_wav_synthesis(complex double *flm, const complex double *f_wav_lm, const complex double *f_scal_lm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min)
{
  int offset, j, l, m;
  int J = s2let_j_max(L, B);
  double wav0, scal0;
  //int l_min = s2let_axisym_el_min(B, J_min);

  offset = 0; 
  for (j = J_min; j <= J; j++){
    for (l = 0; l < L; l++){
      wav0 = sqrt((4.0*PI)/(2.0*l+1.0)) * wav_lm[j*L+l];
      for (m = -l; m <= l; m++){
	flm[lm2ind(l,m)] += f_wav_lm[offset + l*l + l + m] * wav0 ;
      }
    }
    offset += L * L;
  }
  for (l = 0; l < L; l++){
    scal0 = sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lm[l];
    for (m = -l; m <= l; m++){
      flm[lm2ind(l,m)] += f_scal_lm[lm2ind(l,m)] * scal0 ;
    }
  }
}

void s2let_axisym_lm_wav_analysis_multires(complex double *f_wav_lm, complex double *f_scal_lm, const complex double *flm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min)
{
  int bandlimit, offset, j, l, m;
  int J = s2let_j_max(L, B);
  double wav0, scal0;

  offset = 0;
  for (j = J_min; j <= J; j++){
    bandlimit = MIN(s2let_bandlimit(j, J_min, B, L), L);
    for (l = 0; l < bandlimit; l++){
      wav0 = sqrt((4.0*PI)/(2.0*l+1.0)) * wav_lm[j*L+l];
      for (m = -l; m <= l; m++){
	f_wav_lm[offset + l*l + l + m] = flm[lm2ind(l,m)] * wav0 ;
      }
    }
    offset += bandlimit * bandlimit;
  }
  bandlimit = MIN(s2let_bandlimit(J_min-1, J_min, B, L), L);
  for (l = 0; l < bandlimit; l++){
    scal0 = sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lm[l];
    for (m = -l; m <= l; m++){
      f_scal_lm[lm2ind(l,m)] = flm[lm2ind(l,m)] * scal0 ;
    }
  }
}

void s2let_axisym_lm_wav_synthesis_multires(complex double *flm, const complex double *f_wav_lm, const complex double *f_scal_lm, const double *wav_lm, const double *scal_lm, int B, int L, int J_min)
{
  int bandlimit, offset, j, l, m;
  int J = s2let_j_max(L, B);
  double wav0, scal0;
  //int l_min = s2let_axisym_el_min(B, J_min);

  offset = 0; 
  for (j = J_min; j <= J; j++){
    bandlimit = MIN(s2let_bandlimit(j, J_min, B, L), L);
    for (l = 0; l < bandlimit; l++){
      wav0 = sqrt((4.0*PI)/(2.0*l+1.0)) * wav_lm[j*L+l];
      for (m = -l; m <= l; m++){
	flm[lm2ind(l,m)] += f_wav_lm[offset + l*l + l + m] * wav0 ;
      }
    }
    offset += bandlimit * bandlimit;
  }
  bandlimit = MIN(s2let_bandlimit(J_min-1, J_min, B, L), L);
  for (l = 0; l < bandlimit; l++){
    scal0 = sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lm[l];
    for (m = -l; m <= l; m++){
      flm[lm2ind(l,m)] += f_scal_lm[lm2ind(l,m)] * scal0 ;
    }
  }

}

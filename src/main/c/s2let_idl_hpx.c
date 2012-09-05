// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include "idl_export.h"

int s2let_idl_axisym_hpx_wav_analysis_real(int argc, void* argv[])  
{  
  if(argc != 7) return 0;  
  double *f_wav = (double *) argv[0];
  double *f_scal = (double *) argv[1];
  double *f = (double *) argv[2];
  IDL_LONG *nside = (IDL_LONG *) argv[3];
  IDL_LONG *B = (IDL_LONG *) argv[4];
  IDL_LONG *L = (IDL_LONG *) argv[5];
  IDL_LONG *J_min = (IDL_LONG *) argv[6];

  s2let_axisym_hpx_wav_analysis_real(f_wav, f_scal, f, *nside, *B, *L, *J_min);

  return 1;  
}    

int s2let_idl_axisym_hpx_wav_synthesis_real(int argc, void* argv[])  
{  
  if(argc != 7) return 0;  
  double *f = (double *) argv[0];
  double *f_wav = (double *) argv[1];
  double *f_scal = (double *) argv[2];
  IDL_LONG *nside = (IDL_LONG *) argv[3];
  IDL_LONG *B = (IDL_LONG *) argv[4];
  IDL_LONG *L = (IDL_LONG *) argv[5];
  IDL_LONG *J_min = (IDL_LONG *) argv[6];

  s2let_axisym_hpx_wav_synthesis_real(f, f_wav, f_scal, *nside, *B, *L, *J_min);

  return 1;  
}    
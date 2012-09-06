// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"

/*
 * IDL integer types. For historical reasons, we use UCHAR for TYP_BYTE
 * instead of defining an IDL_BYTE type.
 */
 #if defined(ALPHA_OSF) || defined(SUN_64) || defined(LINUX_X86_64) || defined(HPUX_64) || defined(IRIX_64) || defined(AIX_64)
#define IDL_SIZEOF_C_LONG 8
#else
#define IDL_SIZEOF_C_LONG 4
#endif
#if (IDL_SIZEOF_C_LONG == 8) || defined(MSWIN_64)
#define IDL_SIZEOF_C_PTR  8
#else
#define IDL_SIZEOF_C_PTR  4
#endif
typedef short IDL_INT;
typedef unsigned short IDL_UINT;
#if IDL_SIZEOF_C_LONG == 8
typedef int IDL_LONG;
typedef unsigned int IDL_ULONG;
#elif IDL_SIZEOF_C_LONG == 4
typedef long IDL_LONG;
typedef unsigned long IDL_ULONG;
#else
#error "IDL_LONG not defined --- unexpected value of IDL_SIZEOF_C_LONG"
#endif

int s2let_idl_axisym_wav_analysis_real(int argc, void* argv[])  
{  
  if(argc != 6) return 0;  
  double *f_wav = (double *) argv[0];
  double *f_scal = (double *) argv[1];
  double *f = (double *) argv[2];
  IDL_LONG *B = (IDL_LONG *) argv[3];
  IDL_LONG *L = (IDL_LONG *) argv[4];
  IDL_LONG *J_min = (IDL_LONG *) argv[5];

  s2let_axisym_wav_analysis_real(f_wav, f_scal, f, *B, *L, *J_min);

  return 1;  
}    

int s2let_idl_axisym_wav_synthesis_real(int argc, void* argv[])  
{  
  if(argc != 6) return 0;  
  double *f = (double *) argv[0];
  double *f_wav = (double *) argv[1];
  double *f_scal = (double *) argv[2];
  IDL_LONG *B = (IDL_LONG *) argv[3];
  IDL_LONG *L = (IDL_LONG *) argv[4];
  IDL_LONG *J_min = (IDL_LONG *) argv[5];

  s2let_axisym_wav_synthesis_real(f, f_wav, f_scal, *B, *L, *J_min);

  return 1;  
}    

int s2let_idl_axisym_wav_analysis_multires_real(int argc, void* argv[])  
{  
  if(argc != 6) return 0;  
  double *f_wav = (double *) argv[0];
  double *f_scal = (double *) argv[1];
  double *f = (double *) argv[2];
  IDL_LONG *B = (IDL_LONG *) argv[3];
  IDL_LONG *L = (IDL_LONG *) argv[4];
  IDL_LONG *J_min = (IDL_LONG *) argv[5];

  s2let_axisym_wav_analysis_multires_real(f_wav, f_scal, f, *B, *L, *J_min);

  return 1;  
}    

int s2let_idl_axisym_wav_synthesis_multires_real(int argc, void* argv[])  
{  
  if(argc != 6) return 0;  
  double *f = (double *) argv[0];
  double *f_wav = (double *) argv[1];
  double *f_scal = (double *) argv[2];
  IDL_LONG *B = (IDL_LONG *) argv[3];
  IDL_LONG *L = (IDL_LONG *) argv[4];
  IDL_LONG *J_min = (IDL_LONG *) argv[5];

  s2let_axisym_wav_synthesis_multires_real(f, f_wav, f_scal, *B, *L, *J_min);

  return 1;  
}    

// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include <s2let.h>
#include "mex.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))

/**
 * MATLAB interface: s2let_transform_analysis_mw_mex.
 * This function for internal use only.
 * Compute spin directional wavelet transform (analysis)
 * with output in pixel space.
 *
 * Usage:
 *   [f_wav, f_scal] = ...
 *        s2let_transform_analysis_mw_mex(f, B, L, J_min, N, spin, reality, downsample, spin_lowered, original_spin);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int i, j, B, L, J_min, N, spin, f_m, f_n, reality, downsample, normalization, original_spin;
  double *f_wav_real, *f_scal_real, *f_real, *f_wav_imag, *f_scal_imag, *f_imag;
  complex double *f_wav = NULL, *f_scal = NULL, *f = NULL;
  double *f_wav_r = NULL, *f_scal_r = NULL, *f_r = NULL;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=10) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidInput:nrhs",
          "Require ten inputs.");
  }
  if(nlhs!=2) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidOutput:nlhs",
          "Require two outputs.");
  }

  // Parse reality flag
  iin = 6;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidInput:reality",
          "Reality flag must be logical.");
  reality = mxIsLogicalScalarTrue(prhs[iin]);

  // Parse multiresolution flag
  iin = 7;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidInput:downsample",
          "Multiresolution flag must be logical.");
  downsample = mxIsLogicalScalarTrue(prhs[iin]);

  // Parse normalization flag
  iin = 8;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidInput:spinlowered",
          "SpinLowered flag must be logical.");
  if (mxIsLogicalScalarTrue(prhs[iin]))
    normalization = S2LET_WAV_NORM_SPIN_LOWERED;
  else
    normalization = S2LET_WAV_NORM_DEFAULT;

  // Parse original spin
  iin = 9;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidInput:spinloweredfrom",
          "SpinLoweredFrom must be integer.");
  }
  original_spin = (int)mxGetScalar(prhs[iin]);

  // Parse input dataset f
  iin = 0;
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);
  f_real = mxGetPr(prhs[iin]);
  if(reality){
    f_r = (double*)malloc(f_m * f_n * sizeof(double));
    for (i=0; i<f_m*f_n; i++)
      f_r[i] = f_real[i];
  }else{
    f_imag = mxGetPi(prhs[iin]);
    f = (complex double*)malloc(f_m * f_n * sizeof(complex double));
    for (i=0; i<f_m*f_n; i++)
      f[i] = f_real[i] + I * f_imag[i];
  }

  // Parse wavelet parameter B
  iin = 1;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be integer.");
  }
  B = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)B || B <= 1)
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be positive integer greater than 2");

  // Parse harmonic band-limit L
  iin = 2;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit L must be positive integer.");

  if( f_m*f_n != L*(2*L-1) ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidInput:LbandLimit",
          "L must correspond to the sampling scheme, i.e. f = L*(2*L-1) samples.");
  }

  // Parse first scale J_min
  iin = 3;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidInput:Jmin",
          "First scale J_min must be integer.");
  }
  J_min = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)J_min || J_min < 0)
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidInput:Jmin",
          "First scale J_min must be positive integer.");

  // Compute ultimate scale J_max
  int J = s2let_j_max(L, B);

  if( J_min > J+1 ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidInput:Jmin",
          "First scale J_min must be larger than that!");
  }

    // Parse azimuthal/directional band-limit N
  iin = 4;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidInput:NbandLimit",
          "Azimuthal/directional band-limit N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);

    // Parse spin
  iin = 5;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_mw_mex:InvalidInput:spin",
          "spin must be integer.");
  }
  spin = (int)mxGetScalar(prhs[iin]);

  // Perform wavelet transform in harmonic space and then reconstruction.
  if(downsample){
    // Multiresolution algorithm
    if(reality){

    }else{
      s2let_allocate_mw_f_wav_multires(&f_wav, &f_scal, B, L, J_min, N);
      s2let_wav_analysis_mw_multires(f_wav, f_scal, f, B, L, J_min, N, spin, normalization, original_spin);
    }
  }else{
    // Full resolution algorithm
    if(reality){

    }else{
      s2let_allocate_mw_f_wav(&f_wav, &f_scal, B, L, J_min, N);
      s2let_wav_analysis_mw(f_wav, f_scal, f, B, L, J_min, N, spin, normalization, original_spin);
    }
  }

  // Compute size of wavelet array
  int bandlimit, wavsize = 0, scalsize = 0;
  if(downsample){
    for (j = J_min; j <= J; j++){
        bandlimit = MIN(s2let_bandlimit(j, J_min, B, L), L);
        wavsize += (2*N-1) * bandlimit * (2 * bandlimit - 1);
     }
     bandlimit = MIN(s2let_bandlimit(J_min-1, J_min, B, L), L);
     scalsize = bandlimit * (2 * bandlimit - 1);
  }else{
    wavsize = (J+1-J_min) * (2*N-1) * L * ( 2 * L - 1 );
    scalsize = L * ( 2 * L - 1 );
  }

  // Output wavelets
  if(reality){


  }else{

    iout = 0;
    plhs[iout] = mxCreateDoubleMatrix(1, wavsize, mxCOMPLEX);
    f_wav_real = mxGetPr(plhs[iout]);
    f_wav_imag = mxGetPi(plhs[iout]);
    for (i=0; i<wavsize; i++){
      f_wav_real[ i ] = creal( f_wav[ i ] );
      f_wav_imag[ i ] = cimag( f_wav[ i ] );
    }

    iout = 1;
    plhs[iout] = mxCreateDoubleMatrix(1, scalsize, mxCOMPLEX);
    f_scal_real = mxGetPr(plhs[iout]);
    f_scal_imag = mxGetPi(plhs[iout]);
    for (i=0; i<scalsize; i++){
      f_scal_real[i] = creal( f_scal[i] );
      f_scal_imag[i] = cimag( f_scal[i] );
    }

  }

   if(reality){
    //free(f_r);
    //free(f_wav_r);
    //free(f_scal_r);
  }else{
    free(f);
    free(f_wav);
    free(f_scal);
  }

}

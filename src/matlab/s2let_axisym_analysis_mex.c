// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <s2let.h>
#include "mex.h"

/**
 * MATLAB interface: s2let_axisym_analysis.
 * This function for internal use only.
 * Compute axisymmetric wavelet transform (analysis)
 * with output in pixel space.
 *
 * Usage: 
 *   [f_wav, f_scal] = ...
 *        s2let_axisym_analysis_mex(f, B, L, J_min, reality);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int n, i, j, B, L, J_min, f_m, f_n, reality;
  double *f_wav_real, *f_scal_real, *f_real, *f_wav_imag, *f_scal_imag, *f_imag;
  complex double *f_wav = NULL, *f_scal = NULL, *f = NULL;
  double *f_wav_r = NULL, *f_scal_r = NULL, *f_r = NULL;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=5) {
    mexErrMsgIdAndTxt("s2let_axisym_analysis_mex:InvalidInput:nrhs",
          "Require five inputs.");
  }
  if(nlhs!=2) {
    mexErrMsgIdAndTxt("s2let_axisym_analysis_mex:InvalidOutput:nlhs",
          "Require two outputs.");
  }

  // Parse reality flag
  iin = 4;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("s2let_axisym_analysis_mex:InvalidInput:reality",
          "Reality flag must be logical.");
  reality = mxIsLogicalScalarTrue(prhs[iin]);

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
    mexErrMsgIdAndTxt("s2let_axisym_analysis_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be integer.");
  }
  B = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)B || B <= 1)
    mexErrMsgIdAndTxt("s2let_axisym_analysis_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be positive integer greater than 2");

  // Parse harmonic band-limit L
  iin = 2;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_axisym_analysis_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("s2let_axisym_analysis_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit L must be positive integer.");

  if( B >= L ) {
    mexErrMsgIdAndTxt("s2let_axisym_analysis_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be smaller than L!");
  }

  if( f_m*f_n != L*(2*L-1) ) {
    mexErrMsgIdAndTxt("s2let_axisym_analysis_mex:InvalidInput:LbandLimit",
          "L must correspond to the sampling scheme, i.e. f = L*(2*L-1) samples.");
  }
 
  // Parse first scale J_min
  iin = 3;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_axisym_analysis_mex:InvalidInput:Jmin",
          "First scale J_min must be integer.");
  }
  J_min = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)J_min || J_min <= 0)
    mexErrMsgIdAndTxt("s2let_axisym_analysis_mex:InvalidInput:Jmin",
          "First scale J_min must be positive integer.");

  // Compute ultimate scale J_max
  int J = s2let_j_max(L, B);

  if( J_min > J+1 ) {
    mexErrMsgIdAndTxt("s2let_axisym_analysis_mex:InvalidInput:Jmin",
          "First scale J_min must be larger than that!");
  }

  // Perform wavelet transform in harmonic space and then FLAG reconstruction.
  if(reality){
    f_wav_r = (double*)malloc( (J+1) * L * (2*L-1) * sizeof(double));
    f_scal_r = (double*)malloc( L * (2*L-1) * sizeof(double));
    s2let_axisym_wav_analysis_real(f_wav_r, f_scal_r, f_r, B, L, J_min);
  }else{
    f_wav = (complex double*)malloc( (J+1) * L * (2*L-1) * sizeof(complex double));
    f_scal = (complex double*)malloc( L * (2*L-1) * sizeof(complex double));
    s2let_axisym_wav_analysis(f_wav, f_scal, f, B, L, J_min); 
  }


  // Output wavelets
  if(reality){

    iout = 0;
    plhs[iout] = mxCreateDoubleMatrix(J+1, L*(2*L-1), mxREAL);
    f_wav_real = mxGetPr(plhs[iout]);
    for(j=0; j<=J; j++)
      for (i=0; i<L*(2*L-1); i++)
        f_wav_real[ i * (J + 1 ) + j ] = creal(f_wav_r[ j * L*(2*L-1) + i ]);

    iout = 1;
    plhs[iout] = mxCreateDoubleMatrix(1, L*(2*L-1), mxREAL);
    f_scal_real = mxGetPr(plhs[iout]);
    for (i=0; i<f_m*f_n; i++)
      f_scal_real[i] = creal(f_scal_r[i]);

  }else{

    iout = 0;
    plhs[iout] = mxCreateDoubleMatrix(J+1, L*(2*L-1), mxCOMPLEX);
    f_wav_real = mxGetPr(plhs[iout]);
    f_wav_imag = mxGetPi(plhs[iout]);
    for(j=0; j<=J; j++)
      for (i=0; i<L*(2*L-1); i++){
        f_wav_real[ i * (J + 1 ) + j ] = creal( f_wav[ j * L*(2*L-1) + i ] );
        f_wav_imag[ i * (J + 1 ) + j ] = cimag( f_wav[ j * L*(2*L-1) + i ] );
      }

    iout = 1;
    plhs[iout] = mxCreateDoubleMatrix(1, L*(2*L-1), mxCOMPLEX);
    f_scal_real = mxGetPr(plhs[iout]);
    f_scal_imag = mxGetPi(plhs[iout]);
    for (i=0; i<L*(2*L-1); i++){
      f_scal_real[i] = creal( f_scal[i] );
      f_scal_imag[i] = cimag( f_scal[i] );
    }

  }

   if(reality){
    free(f_r);
    free(f_wav_r);
    free(f_scal_r);
  }else{
    free(f);
    free(f_wav);
    free(f_scal);
  }

}

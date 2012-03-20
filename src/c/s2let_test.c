// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <assert.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

double maxerr_cplx(complex double *a, complex double *b, int size)
{
	double value = 0;
	int i;
	for(i = 0; i<size; i++){
		value = MAX( cabs( a[i]-b[i] ), value );
	}
	return value;
}

double maxerr(double *a, double *b, int size)
{
	double value = 0;
	int i;
	for(i = 0; i<size; i++){
		value = MAX( abs( a[i]-b[i] ), value );
	}
	return value;
}

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

void s2let_random_flm(complex double *flm, int L, int seed)
{
	int i;
	srand( time(NULL) );
	for (i=0; i<L*L; i++){
		flm[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
	}
}


void s2let_random_flm_real(complex double *flm, int L, int seed) {

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


void s2let_tilling_test(int B, int L, int J_min)
{		
	double *kappa, *kappa0;
	s2let_allocate_tilling(&kappa, &kappa0, B, L);

	s2let_tilling(kappa, kappa0, B, L, J_min);

	/*
	int l, j;
	int J = s2let_j_max(L, B);
	printf("> KAPPA_0 :       ");
	for (l = 0; l < L; l++){
		printf(" %2.2f,", kappa0[l]);
	}
	printf("\n");
	for (j = 0; j <= J; j++){
		printf("> KAPPA : j = %i : ", j);
		for (l = 0; l < L; l++){
			printf(" %2.2f,", kappa[l+j*L]);
		}
		printf("\n");
	}
	*/

	double res = s2let_check_identity(kappa, kappa0, B, L, J_min);
	printf("  - Identity residuals : %6.5e\n", res);

	free(kappa);
	free(kappa0);
}

void s2let_wav_lm_test(int B, int L, int J_min, int seed)
{
	clock_t time_start, time_end;

	double *wav_lm, *scal_lm;
	s2let_allocate_wav_lm(&wav_lm, &scal_lm, B, L);

	time_start = clock();
	s2let_wav_lm(wav_lm, scal_lm, B, L, J_min);
	time_end = clock();
	printf("  - Generate wavelets  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	/*
	int j, l, m;
	int J = s2let_j_max(L, B);
	printf("-- scal -- \n");
	for (l = 0; l < L; l++){
		printf("scal_lm(%i) = %f\n", l,sqrt((4.0*PI)/(2.0*l+1.0)) * scal_lm[l]);
	}
	for (j = 0; j <= J; j++ ){
		printf("-- j = %i -- \n", j);
		for (l = 0; l < L; l++){
			printf("wav_lm(%i) = %f\n", l, sqrt((4.0*PI)/(2.0*l+1.0)) * wav_lm[j*L+l]);
		}
	}
	*/	
	 
	complex double *f_wav_lm, *f_scal_lm, *flm, *flm_rec;
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	flm_rec = (complex double*)calloc(L * L, sizeof(complex double));

	s2let_random_flm(flm, L, seed);
	
	s2let_allocate_f_wav_lm(&f_wav_lm, &f_scal_lm, B, L);
	
	time_start = clock();
	s2let_wav_analysis_lm(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);


	time_start = clock();
	s2let_wav_synthesis_lm(flm_rec, f_wav_lm, f_scal_lm, wav_lm, scal_lm, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	/*
	for (l = 0; l < L; l++){
	    for (m = -l; m <= l ; m++){
		    int ind = l*l+l+m;
			printf("(l,m) = (%i,%i) - flm = (%f,%f) - rec = (%f,%f)\n",l,m,creal(flm[ind]),cimag(flm[ind]),creal(flm_rec[ind]),cimag(flm_rec[ind]));
	    }
	}
	*/

	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flm, flm_rec, L*L));

	free(flm);
	free(flm_rec);
	free(f_wav_lm);
	free(f_scal_lm);
	free(wav_lm);
	free(scal_lm);
}

void s2let_wav_test(int B, int L, int J_min, int seed)
{
	clock_t time_start, time_end;
	int spin = 0;
	int verbosity = 0;
	ssht_dl_method_t dl_method = SSHT_DL_RISBO;
	int J = s2let_j_max(L, B);

	complex double *f, *f_rec, *flm, *flm_rec;
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	flm_rec = (complex double*)calloc(L * L, sizeof(complex double));
	f = (complex double*)calloc(L*(2*L-1), sizeof(complex double));
	f_rec = (complex double*)calloc(L*(2*L-1), sizeof(complex double));

	s2let_random_flm(flm, L, seed);

	ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

	complex double *f_wav, *f_scal;
	f_scal = (complex double*)calloc(L*(2*L-1), sizeof(complex double));
	f_wav = (complex double*)calloc((J+1)*L*(2*L-1), sizeof(complex double));

	time_start = clock();
	s2let_wav_analysis(f_wav, f_scal, f, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	time_start = clock();
	s2let_wav_synthesis(f_rec, f_wav, f_scal, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	ssht_core_mw_forward_sov_conv_sym(flm_rec, f_rec, L, spin, dl_method, verbosity);

	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flm, flm_rec, L*L));

	/*
	int l, m;
	for (l = 0; l < L; l++){
	    for (m = -l; m <= l ; m++){
		    int ind = l*l+l+m;
			printf("(l,m) = (%i,%i) - flm = (%f,%f) - rec = (%f,%f)\n",l,m,creal(flm[ind]),cimag(flm[ind]),creal(flm_rec[ind]),cimag(flm_rec[ind]));
	    }
	}
	*/
	
	free(f);
	free(f_rec);
	free(f_wav);
	free(f_scal);
}

void s2let_wav_real_test(int B, int L, int J_min, int seed)
{
	clock_t time_start, time_end;
	int verbosity = 0;
	ssht_dl_method_t dl_method = SSHT_DL_RISBO;
	int J = s2let_j_max(L, B);

	complex *flm, *flm_rec;
	double *f, *f_rec;
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	flm_rec = (complex double*)calloc(L * L, sizeof(complex double));
	f = (double*)calloc(L*(2*L-1), sizeof(double));
	f_rec = (double*)calloc(L*(2*L-1), sizeof(double));

	s2let_random_flm_real(flm, L, seed);

	ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

	double *f_wav, *f_scal;
	f_scal = (double*)calloc(L*(2*L-1), sizeof(double));
	f_wav = (double*)calloc((J+1)*L*(2*L-1), sizeof(double));

	time_start = clock();
	s2let_wav_analysis_real(f_wav, f_scal, f, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet analysis   : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	time_start = clock();
	s2let_wav_synthesis_real(f_rec, f_wav, f_scal, B, L, J_min);
	time_end = clock();
	printf("  - Wavelet synthesis  : %4.4f seconds\n", 
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	ssht_core_mw_forward_sov_conv_sym_real(flm_rec, f_rec, L, dl_method, verbosity);

	printf("  - Maximum abs error  : %6.5e\n", 
		maxerr_cplx(flm, flm_rec, L*L));

	/*
	int l, m;
	for (l = 0; l < L; l++){
	    for (m = -l; m <= l ; m++){
		    int ind = l*l+l+m;
			printf("(l,m) = (%i,%i) - flm = (%f,%f) - rec = (%f,%f)\n",l,m,creal(flm[ind]),cimag(flm[ind]),creal(flm_rec[ind]),cimag(flm_rec[ind]));
	    }
	}
	*/
	
	free(f);
	free(f_rec);
	free(f_wav);
	free(f_scal);
}

int main(int argc, char *argv[]) 
{
	
	const int L = 64;
	const int B = 3;
	const int J_min = 3;
	const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);
	int l_min = s2let_el_min(B, J_min);

	printf("==========================================================\n");
	printf("PARAMETERS: ");
	printf("  L = %i   B = %i   l_wav_min = %i   seed = %i\n", L, B, l_min, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing harmonic tilling...\n");
	s2let_tilling_test(B, L, J_min);
	printf("----------------------------------------------------------\n");
	printf("> Testing axisymmetric wavelets in harmonics space...\n");
	s2let_wav_lm_test(B, L, J_min, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing axisymmetric wavelets in pixel space...\n");
	s2let_wav_test(B, L, J_min, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing real axisymmetric wavelets in pixel space...\n");
	s2let_wav_real_test(B, L, J_min, seed);
	printf("==========================================================\n");


	return 0;		
}

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

void s2let_tilling_test(int B, int L)
{
	int l, j;
	int J = s2let_j_max(L, B);
		
	double *kappa, *kappa0;
	allocate_tilling(&kappa, &kappa0, B, L);

	s2let_tilling(kappa, kappa0, B, L);

	//printf("\n");
	j = 0;
	printf("> KAPPA_0 :       ");
	for (l = 0; l < L; l++){
		printf(" %2.2f,", kappa0[l+j*L]);
	}
	printf("\n");
	for (j = 0; j <= J; j++){
		printf("> KAPPA : j = %i : ", j);
		for (l = 0; l < L; l++){
			printf(" %2.2f,", kappa[l+j*L]);
		}
		printf("\n");
	}

	s2let_check_identity(kappa, kappa0, B, L);
	printf("\n");
}

void s2let_wav_lm_test(int B, int L, int seed)
{
	int j, l, m;
	int J = s2let_j_max(L, B);

	double *wav_lm, *scal_lm;
	s2let_allocate_wav_lm_real(&wav_lm, &scal_lm, B, L);

	int J_min = 0;
	s2let_wav_lm(wav_lm, scal_lm, B, L);

	
	printf("-- scal -- \n");
	for (l = 0; l < L; l++){
		for (m = -l; m <= l; m++){
			printf("scal_lm(%i,%i) = %f\n", l, m,scal_lm[lm2ind(l,m)]);
		}
	}
	for (j = 0; j <= J; j++ ){
		printf("-- j = %i -- \n", j);
		for (l = 0; l < L; l++){
			for (m = -l; m <= l; m++){
				printf("wav_lm(%i,%i) = %f\n", l, m, wav_lm[jlm2ind(j,l,m,L)]);
			}
		}
	}
	 

	complex double *f_wav_lm, *f_scal_lm, *flm, *flm_rec;
	flm = (complex double*)calloc(L * L, sizeof(complex double));
	flm_rec = (complex double*)calloc(L * L, sizeof(complex double));
	s2let_random_flm(flm, L, seed);
	printf("Allocate f_wav\n");
	s2let_allocate_wav_lm(&f_wav_lm, &f_scal_lm, B, L);
	printf("Analysis\n");
	s2let_wav_analysis_lm(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, B, L, J_min);
	printf("Synthesis\n");
	s2let_wav_synthesis_lm(flm_rec, f_wav_lm, f_scal_lm, wav_lm, scal_lm, B, L, J_min);
	for (l = 0; l < L; l++){
	for (m = -l; m <= l ; m++){
		 int ind = l*l+l+m;
			printf("(l,m) = (%i,%i) - flm = (%f,%f) - rec = (%f,%f)\n",l,m,creal(flm[ind]),cimag(flm[ind]),creal(flm_rec[ind]),cimag(flm_rec[ind]));
	}}
	printf("  - Maximum absolute error    : %6.5e\n", 
		maxerr_cplx(flm, flm_rec, L*L));

}

int main(int argc, char *argv[]) 
{
	
	const int L = 12;
	const int B = 2;
	const int seed = (int)(1000000.0*(double)clock()/(double)CLOCKS_PER_SEC)/10;

	printf("=========================================================\n");
	printf("PARAMETERS : ");
	printf("  L = %i   B = %i  \n", L, B);
	printf("----------------------------------------------------------\n");
	printf("> Testing harmonic tilling...\n");
	s2let_tilling_test(B, L);
	printf("OK\n");
	printf("----------------------------------------------------------\n");
	printf("> Testing wavelets in harmonics space...\n");
	s2let_wav_lm_test(B, L, seed);
	printf("OK\n");

	return 0;		
}

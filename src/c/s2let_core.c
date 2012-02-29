// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))

// Computes minimum harmonic index supported by needlets.
int s2let_el_min(int B, int J_min)
{
	return ceil(pow(B, J_min));
}

// Computes needlet maximum level required to ensure exact reconstruction.
int s2let_j_max(int L, int B)
{
	return ceil(log(L) / log(B));
}

// Computes needlet band-limits for all levels j.
void s2let_Lj(int *Lj, int B, int J_min, int L)
{
	int j;
	int J = s2let_j_max(L, B);
	for (j = J_min ; j <= J ; j++){
		Lj[j+1-J_min] = MIN(L, ceil(pow(B,j+1)) + 1);
	}
}

void allocate_tilling(double **kappa, double **kappa0, int B, int L)
{
	int J = s2let_j_max(L, B);
	*kappa = (double*)calloc((J+1) * L, sizeof(double));
	*kappa0 = (double*)calloc(L, sizeof(double));
}

void s2let_tilling(double *kappa, double *kappa0, int B, int L)
{
	int j, l;
	int J = s2let_j_max(L, B);

	double *phi2 = (double*)calloc((J+2)*L, sizeof(double));

	int n = 50;

	double kappanorm = s2let_kappa0_quadtrap(1.0 / (double)B, 1.0, n, B);
	//printf("\nKAPPANORM = %f\n", kappanorm);

	for (j = 0; j <= J+1; j++){
		for (l = 0; l < L; l++){
			if (l < pow(B,j-1)) {
				phi2[l+j*L] = 1;
			} else if (l > pow(B,j)) {
				phi2[l+j*L] = 0;
			} else {
				phi2[l+j*L] = s2let_kappa0_quadtrap((double)l/pow(B, j), 1.0, n, B) / kappanorm;
			}
		}
	}

	for (l = 0; l < L; l++){
		kappa0[l] = sqrt(phi2[l]);
	}
	for (j = 0; j <= J; j++){
		for (l = 0; l < L; l++){
			kappa[l+j*L] = sqrt(phi2[l+(j+1)*L] - phi2[l+j*L]);
		}
	}

}


void s2let_check_identity(double *kappa, double *kappa0, int B, int L)
{
	int l, j;
	int J = s2let_j_max(L, B);
		
	double *ident;
	ident = (double*)calloc(L, sizeof(double));

	printf("Tilling identity: ");
	ident[0] = pow(kappa0[0], 2.0);
	printf(" %2.2f ", ident[0]);
	for (l = 1; l < L; l++){
		for (j = 0; j <= J; j++){
			ident[l] += pow(kappa[l+j*L], 2.0);
		}
		printf(" %2.2f ", ident[l]);
	}

}

int jlm2ind(int j, int l, int m, int L)
{
	return j*L*L + l*l + l + m ;
}

void s2let_allocate_wav_lm(double **wav_lm, double **scal_lm, int B, int L)
{
	int J = s2let_j_max(L, B);
	*wav_lm = (double*)calloc((J+1) * L, sizeof(double));
	*scal_lm = (double*)calloc(L, sizeof(double));
}

void s2let_wav_lm(double *wav_lm, double *scal_lm, int B, int L, int J_min)
{
	int j, l, m;
	int J = s2let_j_max(L, B);
	int l_min = s2let_el_min(B, J_min);

	double *kappa, *kappa0;
	allocate_tilling(&kappa, &kappa0, B, L);
	s2let_tilling(kappa, kappa0, B, L);

	for (j = J_min; j <= J; j++){
		for (l = l_min; l < L; l++){
			wav_lm[j * L + l] = sqrt( (2 * l + 1) / (4.0 * PI) ) * kappa[l + j * L];
		}
	}
	for (l = 0; l < l_min; l++){
		scal_lm[l] = sqrt( (2 * l + 1) / (4.0 * PI) ) * kappa0[l];
	}

}

void s2let_allocate_f_wav_lm(double **f_wav_lm, double **f_scal_lm, int B, int L)
{
	int J = s2let_j_max(L, B);
	*f_wav_lm = (double*)calloc((J+1) * L * L, sizeof(double));
	*f_scal_lm = (double*)calloc(L * L, sizeof(double));
}

void s2let_wav_analysis_lm(complex double *f_wav_lm, complex double *f_scal_lm, const complex double *flm, double *wav_lm, double *scal_lm, int B, int L, int J_min)
{
	int j, l, m;
	int J = s2let_j_max(L, B);
	double wav0, scal0;
	int l_min = s2let_el_min(B, J_min);

	for (j = J_min; j <= J; j++){
		for (l = l_min; l < L; l++){
			wav0 = wav_lm[j * L + l];
			for (m = -l; m <= l; m++){
				f_wav_lm[jlm2ind(j,l,m,L)] = wav0;
			}
		}
	}
	for (l = 0; l < l_min; l++){
		scal0 = scal_lm[l];
		for (m = -l; m <= l; m++){
			f_scal_lm[l*l + l + m] = scal0;
		}
	}
}

void s2let_wav_analysis_lm_real(double *f_wav_lm, double *f_scal_lm, const double *flm, double *wav_lm, double *scal_lm, int B, int L, int J_min)
{	
	int j, l, m;
	int J = s2let_j_max(L, B);
	double wav0, scal0;
	int l_min = s2let_el_min(B, J_min);

	for (j = J_min; j <= J; j++){
		for (l = l_min; l < L; l++){
			wav0 = wav_lm[j * L + l];
			for (m = -l; m <= l; m++){
				f_wav_lm[jlm2ind(j,l,m,L)] = wav0;
			}
		}
	}
	for (l = 0; l < l_min; l++){
		scal0 = scal_lm[l];
		for (m = -l; m <= l; m++){
			f_scal_lm[l*l + l + m] = scal0;
		}
	}	
}

void s2let_allocate_f_wav(double **f_wav_lm, double **f_scal_lm, int B, int L)
{
	int J = s2let_j_max(L, B);
	*f_wav_lm = (double*)calloc((J+1) * L * L, sizeof(double));
	*f_scal_lm = (double*)calloc(L * L, sizeof(double));
}

void s2let_wav_analysis(complex double *f_wav, complex double *f_scal, const complex double *f, int B, int L, int J_min){
	
}

void s2let_wav_analysis_real(double *f_wav, double *f_scal, const complex double *f, int B, int L, int J_min){
	
}





// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <math.h>
#include <stdlib.h>

int s2let_bandlimit(int B, int j)
{
	return ceil(pow(B, j+1));
}

int s2let_el_min(int B, int J_min)
{
	return ceil(pow(B, J_min));
}

int s2let_j_max(int L, int B)
{
	return ceil(log(L) / log(B));
}

void s2let_axisym_allocate_tiling(double **kappa, double **kappa0, int B, int L)
{
	int J = s2let_j_max(L, B);
	*kappa = (double*)calloc((J+1) * L, sizeof(double));
	*kappa0 = (double*)calloc(L, sizeof(double));
}

void s2let_axisym_tiling(double *kappa, double *kappa0, int B, int L, int J_min)
{
	int j, l;
	int J = s2let_j_max(L, B);

	double previoustemp, temp;
	double *phi2 = (double*)calloc((J+2) * L, sizeof(double));

	s2let_tiling_phi2(phi2, B, L, J_min);

	for (l = 0; l < L; l++){
		kappa0[l] = sqrt(phi2[l+J_min*L]);
	}
	
	for (j = J_min; j <= J; j++){
		for (l = 0; l < L; l++){
			temp = sqrt(phi2[l+(j+1)*L] - phi2[l+j*L]);
			if( isnan(temp) || isinf(temp) )
				kappa[l+j*L] = previoustemp;
			else
				kappa[l+j*L] = temp;
			previoustemp = temp;
		}
		for (l = 0; l < L; l++){
			if( !finite(kappa[l+j*L]) ) 
				 kappa[l+j*L] = kappa[l+j*L-1];
		}
	}

	free(phi2);
}

void s2let_tiling_phi2(double *phi2, int B, int L, int J_min)
{
	int j, l;
	int J = s2let_j_max(L, B);
	int n = 300;

	double kappanorm = s2let_kappa0_quadtrap(1.0 / (double)B, 1.0, n, B);
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

}

double s2let_axisym_check_identity(double *kappa, double *kappa0, int B, int L, int J_min)
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

// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "s2let.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))

/*!
 * Computes minimum harmonic index supported by needlets.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval ell_min
 */
int s2let_el_min(int B, int J_min)
{
	return ceil(pow(B, J_min));
}

/*!
 * Computes needlet maximum level required to ensure exact reconstruction.
 *
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  B Wavelet parameter.
 * \retval j_max
 */
int s2let_j_max(int L, int B)
{
	return ceil(log(L) / log(B));
}

// Computes needlet band-limits for all levels j.
void s2let_Lj(int *Lj, int B, int L, int J_min)
{
	int j;
	int J = s2let_j_max(L, B);
	for (j = J_min ; j <= J ; j++){
		Lj[j+1-J_min] = MIN(L, ceil(pow(B,j+1)) + 1);
	}
}

void s2let_allocate_tilling(double **kappa, double **kappa0, int B, int L)
{
	int J = s2let_j_max(L, B);
	*kappa = (double*)calloc((J+1) * L, sizeof(double));
	*kappa0 = (double*)calloc(L, sizeof(double));
}

/*!
 * Generates tilling in harmonic space.
 *
 * \param[out]  kappa Wavelet kernel functions for wavelets.
 * \param[out]  kappa0 Kernel for the scaling function.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_tilling(double *kappa, double *kappa0, int B, int L, int J_min)
{
	int j, l;
	int J = s2let_j_max(L, B);

	double *phi2 = (double*)calloc((J+2) * L, sizeof(double));

	s2let_tilling_phi2(phi2, B, L, J_min);

	for (l = 0; l < L; l++){
		kappa0[l] = sqrt(phi2[l+J_min*L]);
	}
	
	for (j = J_min; j <= J; j++){
		for (l = 0; l < L; l++){
			kappa[l+j*L] = sqrt(phi2[l+(j+1)*L] - phi2[l+j*L]);
		}
	}

	free(phi2);
}

/*!
 * Generates smooth functions to construct the tilling.
 *
 * \param[out]  phi2 Smooth step functions for wavelets.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_tilling_phi2(double *phi2, int B, int L, int J_min)
{
	int j, l;
	int J = s2let_j_max(L, B);

	int n = 100;

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

/*!
 * Checks exactness of the harmonic tilling.
 *
 * \param[in]  kappa Wavelet kernel functions for wavelets.
 * \param[in]  kappa0 Kernel for the scaling function.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval Achieved accuracy (should be 0).
 */
double s2let_check_identity(double *kappa, double *kappa0, int B, int L, int J_min)
{
	int l, j;
	int J = s2let_j_max(L, B);
	int l_min = s2let_el_min(B, J_min);
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

/*!
 * Indice corresponding to a triplet (j, l, m) in the wavelets.
 *
 * \param[in]  j Scale indice.
 * \param[in]  l Multipole indice.
 * \param[in]  m Order indice.
 * \param[in]  L Angular harmonic band-limit.
 * \retval Indice
 */
int jlm2ind(int j, int l, int m, int L)
{
	return j*L*L + l*l + l + m ;
}

/*!
 * Indice corresponding to a doublet (l, m).
 *
 * \param[in]  l Multipole indice.
 * \param[in]  m Order indice.
 * \param[in]  L Angular harmonic band-limit.
 * \retval Indice
 */
int lm2ind(int l, int m)
{
	return l*l + l + m ;
}




# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# if you want to use the Numpy-C-API from Cython
# (not strictly necessary for this example)
np.import_array()

cdef extern from "s2let.h":

	void s2let_transform_axisym_lm_allocate_wav(
		double **wav_lm, double **scal_lm, const s2let_parameters_t *parameters);

	void s2let_transform_axisym_lm_wav(
		double *wav_lm, double *scal_lm, const s2let_parameters_t *parameters);

	void s2let_transform_axisym_lm_wav_analysis(
		double complex *f_wav_lm,
		double complex *f_scal_lm,
		const double complex *flm,
		const double *wav_lm,
		const double *scal_lm,
		const s2let_parameters_t *parameters
	);

	void s2let_transform_axisym_lm_wav_synthesis(
		double complex *flm,
		const double complex *f_wav_lm,
		const double complex *f_scal_lm,
		const double *wav_lm,
		const double *scal_lm,
		const s2let_parameters_t *parameters
	);

	ctypedef struct s2let_parameters_t:
		int B
		int L
		int J_min

	int s2let_j_max(s2let_parameters_t *parameters)

cdef extern from "stdlib.h":
	void free(void* ptr)

def pys2let_j_max(B, L, J_min):
	cdef s2let_parameters_t parameters = {};
	parameters.B = B;
	parameters.L = L;
	parameters.J_min = J_min;
	return s2let_j_max(&parameters);

def healpy_lm(el, em, L):
	return em*(2*L-1-em)/2+el

def pys2let_transform_axisym_lm_wav_analysis(
	np.ndarray[double complex, ndim=1, mode="c"] flm_hp not None, B, L, J_min):

	cdef s2let_parameters_t parameters = {};
	parameters.B = B;
	parameters.L = L;
	parameters.J_min = J_min;
	J = s2let_j_max(&parameters);

	cdef double *wav_lm, *scal_lm;
	s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, &parameters);
	s2let_transform_axisym_lm_wav(wav_lm, scal_lm, &parameters);

	f_lm = np.empty([L * L,], dtype=complex)
	f_scal_lm = np.empty([L * L,], dtype=complex)
	f_wav_lm = np.empty([L * L * (J+1-J_min),], dtype=complex)

	for el from 0 <= el < L:
		for em from 0 <= em <= el:
			f_lm[ el * el + el - em ] = pow(-1.0, -em) * ( flm_hp[healpy_lm(el, em, L)] ).conjugate()
			f_lm[ el * el + el + em ] = flm_hp[healpy_lm(el, em, L)]

	s2let_transform_axisym_lm_wav_analysis(
		<double complex*> np.PyArray_DATA(f_wav_lm),
		<double complex*> np.PyArray_DATA(f_scal_lm),
		<double complex*> np.PyArray_DATA(f_lm),
		wav_lm, scal_lm, &parameters
	);

	free(scal_lm);
	free(wav_lm);

	f_scal_lm_hp = np.empty([L*(L+1)/2,], dtype=complex)
	f_wav_lm_hp = np.empty([L*(L+1)/2, J-J_min+1], dtype=complex)

	for el from 0 <= el < L:
		for em from 0 <= em <= el:
			f_scal_lm_hp[healpy_lm(el, em, L)] = f_scal_lm[ el * el + el + em ]
			for j from 0 <= j <= J-J_min:
				f_wav_lm_hp[healpy_lm(el, em, L), j] = f_wav_lm[ L*L*j + el * el + el + em ]

	return f_wav_lm_hp, f_scal_lm_hp


def pys2let_transform_axisym_lm_wav_synthesis(
	np.ndarray[double complex, ndim=2, mode="c"] f_wav_lm_hp not None, 
	np.ndarray[double complex, ndim=1, mode="c"] f_scal_lm_hp not None, B, L, J_min):

	cdef s2let_parameters_t parameters = {};
	parameters.B = B;
	parameters.L = L;
	parameters.J_min = J_min;
	J = s2let_j_max(&parameters);

	cdef double *wav_lm, *scal_lm;
	s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, &parameters);
	s2let_transform_axisym_lm_wav(wav_lm, scal_lm, &parameters);

	f_lm = np.empty([L * L,], dtype=complex)
	f_scal_lm = np.empty([L * L,], dtype=complex)
	f_wav_lm = np.empty([L * L * (J+1-J_min),], dtype=complex)

	for el from 0 <= el < L:
		for em from 0 <= em <= el:
			f_scal_lm[ el * el + el - em ] = pow(-1.0, -em) * f_scal_lm_hp[healpy_lm(el, em, L)].conjugate()
			f_scal_lm[ el * el + el + em ] = f_scal_lm_hp[healpy_lm(el, em, L)]
			for j from 0 <= j <= J-J_min:
				f_wav_lm[ L*L*j + el * el + el - em ] = pow(-1.0, -em) * f_wav_lm_hp[healpy_lm(el, em, L), j].conjugate()
				f_wav_lm[ L*L*j + el * el + el + em ] = f_wav_lm_hp[healpy_lm(el, em, L), j]


	s2let_transform_axisym_lm_wav_synthesis(
		<double complex*> np.PyArray_DATA(f_lm),
		<double complex*> np.PyArray_DATA(f_wav_lm),
		<double complex*> np.PyArray_DATA(f_scal_lm),
		wav_lm, scal_lm, &parameters
	);

	free(scal_lm);
	free(wav_lm);

	f_lm_hp = np.empty([L*(L+1)/2,], dtype=complex)

	for el from 0 <= el < L:
		for em from 0 <= em <= el:
			f_lm_hp[healpy_lm(el, em, L)] = f_lm[ el * el + el + em ]

	return f_lm_hp

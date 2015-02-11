
import healpy as hp
from pys2let import *
import math
import matplotlib.pyplot as plt

nside = 128
L = 128
J_min = 1
B = 3
J = pys2let_j_max(B, L, J_min)

# The filename of some random healpix map
fname = '/Users/bl/Dropbox/Astrodata/SDSS/Fields/Planck_EBV_256rQ.fits'
#fname = '/Users/bl/Dropbox/Wavelets/s2let/data/somecmbsimu_hpx_128.fits'

# Read healpix map and compute alms. 
# f_lm has size L*(L+1)/2
f = hp.read_map(fname)
f_lm = hp.map2alm(f, lmax=L-1)

# Call pys2let and compute wavelet transform. Returns the harmonic coefficients of the wavelets.
# f_scal_lm has size L*(L+1)/2
# f_wav_lm has size L*(L+1)/2 by J-J_min+1
f_wav_lm, f_scal_lm = pys2let_transform_axisym_lm_wav_analysis(f_lm, B, L, J_min)

# Reconstruct healpix maps on the sphere and plot them
f_scal = hp.alm2map(f_scal_lm, nside=nside, lmax=L-1)
hp.mollview(f_scal)
f_wav = np.empty([12*nside*nside, J-J_min+1])
for j in range(J-J_min+1):
	flm = f_wav_lm[:,j].ravel()
	f_wav[:,j] = hp.alm2map(flm, nside=nside, lmax=L-1)
	hp.mollview(f_wav[:,j])

plt.show()
% s2let_hpxtest
% Run test for HEALPIX data.
% The wavelet transforms will not be exact but should
% at least reconstruct the input maps with correct precision.
%
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

clear all;
close all;

setenv('HEALPIX','/Users/bl/software/Healpix_2.20a')
setenv('HEALPIXDATA','/Users/bl/software/Healpix_2.20a/data')

% Main parameters
L = 128;
nside = 128;
B = 2;
J_min = 1;
J = s2let_jmax(L, B);
npix = 12*nside*nside;


disp('Generates band-limited function')
flm = zeros(L^2,1); 
flm = rand(size(flm)) + sqrt(-1)*rand(size(flm));
flm = 2.*(flm - (1+sqrt(-1))./2);

disp('Constraint on flms to generate real signal')
for el = 0:L-1
   ind = el*el + el + 1;
   flm(ind,1) = real(flm(ind,1));
   for m = 1:el
      ind_pm = el*el + el + m + 1;
      ind_nm = el*el + el - m + 1;
      flm(ind_nm,1) = (-1)^m * conj(flm(ind_pm,1));
   end  
end


disp('Perform spherical harmonic decomposition with default parameters')
f = s2let_hpx_alm2map(flm, nside, 'L', L);
flm_rec = s2let_hpx_map2alm(f, 'L', L);
default = max(max(abs(flm-flm_rec)))


disp('Perform spherical harmonic decomposition with custom parameters')
f = s2let_hpx_alm2map(flm, nside, 'L', L);
flm_rec = s2let_hpx_map2alm(f, 'nside', nside, 'L', L);
custom = max(max(abs(flm-flm_rec)))
 

disp('Perform axisym wavelet transform with default parameters')
[f_wav, f_scal] = s2let_axisym_hpx_analysis(f);
f_rec = s2let_axisym_hpx_synthesis(f_wav, f_scal);
default = max(max(abs(f-f_rec)))


disp('Perform axisym wavelet transform with custom parameters')
[f_wav, f_scal] = s2let_axisym_hpx_analysis(f, 'B', B, 'L', L, 'J_min', J_min);
f_rec = s2let_axisym_hpx_synthesis(f_wav, f_scal, 'B', B, 'L', L, 'J_min', J_min);
custom = max(max(abs(f-f_rec)))
% s2let_fulltest - Run all tests
%
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

clear all;
close all;

% Main parameters
L = 128;
B = 3;
J_min = 4;
J = s2let_jmax(L, B)

% Checks tilling of harmonic space for axysimmetric wavelets
[kappa kappa0] = s2let_axisym_tilling(B, L, J_min);
error_on_axisym_tilling = s2let_check_axisym_tilling(kappa, kappa0, L, J)

% Generates band-limited function
flm = zeros(L^2,1); 
flm = rand(size(flm)) + sqrt(-1)*rand(size(flm));
flm = 2.*(flm - (1+sqrt(-1))./2);
% Real space synthesis
f = ssht_inverse(flm, L, 'Method', 'MW');

% Perform axisym wavelet transform with default parameters
[f_wav, f_scal] = s2let_axisym_analysis(f);
f_rec = s2let_axisym_synthesis(f_wav, f_scal);
default = max(max(abs(f-f_rec)))

% Perform axisym wavelet transform with default parameters
[f_wav, f_scal] = s2let_axisym_analysis(f, 'downsample', true);
f_rec = s2let_axisym_synthesis(f_wav, f_scal, 'downsample', true);
default_multires = max(max(abs(f-f_rec)))

% Perform axisym wavelet transform with default parameters
[f_wav, f_scal] = s2let_axisym_analysis(f, 'downsample', false);
f_rec = s2let_axisym_synthesis(f_wav, f_scal, 'downsample', false);
default_fullres = max(max(abs(f-f_rec)))

% Perform axisym wavelet transform 
[f_wav, f_scal] = s2let_axisym_analysis(f, 'B', B, 'L', L, 'J_min', J_min);
f_rec = s2let_axisym_synthesis(f_wav, f_scal, 'B', B, 'L', L, 'J_min', J_min);
custom = max(max(abs(f-f_rec)))

% Constraint on flm's to generate real signal
for el = 0:L-1
   ind = el*el + el + 1;
   flm(ind,1) = real(flm(ind,1));
   for m = 1:el
      ind_pm = el*el + el + m + 1;
      ind_nm = el*el + el - m + 1;
      flm(ind_nm,1) = (-1)^m * conj(flm(ind_pm,1));
   end  
end
% Real space synthesis
f_real = ssht_inverse(flm, L, 'Method', 'MW', 'Reality', true);


% Perform real axisym wavelet transform 
[f_wav_real, f_scal_real] = s2let_axisym_analysis(f_real, 'Reality', true);
f_real_rec = s2let_axisym_synthesis(f_wav_real, f_scal_real, 'Reality', true);
default = max(max(abs(f_real-f_real_rec)))

% Perform real axisym wavelet transform 
[f_wav_real, f_scal_real] = s2let_axisym_analysis(f_real, 'downsample', true, 'Reality', true);
f_real_rec = s2let_axisym_synthesis(f_wav_real, f_scal_real, 'downsample', true, 'Reality', true);
default_multires = max(max(abs(f_real-f_real_rec)))

% Perform real axisym wavelet transform 
[f_wav_real, f_scal_real] = s2let_axisym_analysis(f_real, 'downsample', false, 'Reality', true);
f_real_rec = s2let_axisym_synthesis(f_wav_real, f_scal_real, 'downsample', false, 'Reality', true);
default_fullres = max(max(abs(f_real-f_real_rec)))

% Perform real axisym wavelet transform 
[f_wav_real, f_scal_real] = s2let_axisym_analysis(f_real, 'B', B, 'L', L, 'J_min', J_min, 'Reality', true);
f_real_rec = s2let_axisym_synthesis(f_wav_real, f_scal_real, 'B', B, 'L', L, 'J_min', J_min, 'Reality', true);
custom = max(max(abs(f_real-f_real_rec)))


% s2let_fulltest - Run all tests
%
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

clear all;
close all;

% Main parameters
L = 32
B = 3
J_min = 2
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

% Perform axisym wavelet transform 
[f_wav, f_scal] = s2let_axisym_analysis(f, B, L, J_min);
f_rec = s2let_axisym_synthesis(f_wav, f_scal, B, L, J_min);
s2let_transform_error = max(max(abs(f-f_rec)))

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
[f_wav_real, f_scal_real] = s2let_axisym_analysis(f_real, B, L, J_min, 'Reality', true);
f_real_rec = s2let_axisym_synthesis(f_wav_real, f_scal_real, B, L, J_min, 'Reality', true);
s2let_real_transform_error = max(max(abs(f_real-f_real_rec)))


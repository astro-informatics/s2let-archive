% s2let_fulltest - Run all tests
%
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

clear all;
close all;

% Main parameters
L = 32
B = 2
J_min = 2
J = s2let_jmax(L, B)

% Checks tilling of harmonic space for axysimmetric wavelets
[kappa kappa0] = s2let_axisym_tilling(B, L, J_min);
s2let_check_axisym_tilling(kappa, kappa0, L, J);
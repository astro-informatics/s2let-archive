function [kappa kappa0] = s2let_axisym_tilling(B, L, J_min)

% s2let_axisym_tilling - Compute tilling in harmonic space.
% -- Axisymmetric wavelets on the sphere.
%
% Default usage :
%
%   [kappa kappa0] = s2let_axisym_tilling(B, L, J_min)
%
% kappa is an array containing wavelet tilling .
% kappa0 contains the scaling function.
% B is the wavelet parameter,
% L is the angular band-limit,
% J_min the first wavelet to be used.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

p = inputParser;
p.addRequired('B', @isnumeric);          
p.addRequired('L', @isnumeric);   
p.addRequired('J_min', @isnumeric); 
p.parse(B, L, J_min);
args = p.Results;

[kappa kappa0] = s2let_axisym_tilling_mex(B, L, J_min);
  
end
function [psi_lm phi_l] = s2let_wavelet_tiling(B, L, N, Spin, J_min)

% s2let_wavelet_tiling - Compute tiling in harmonic space.
% -- Spin directional wavelets on the sphere.
%
% Default usage :
%
%   [psi_lm phi_l] = s2let_wavelet_tiling(B, L, N, Spin, J_min)
%
% psi_lm is an array containing the wavelet spherical harmonic coefficients.
% phi_l is an array containing the scaling function spherical harmonic coefficients (l only).
% B is the wavelet parameter,
% L is the angular band-limit,
% N is the azimuthal/directional band-limit,
% Spin is the spin number,
% J_min the first wavelet to be used.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

p = inputParser;
p.addRequired('B', @isnumeric);
p.addRequired('L', @isnumeric);
p.addRequired('N', @isnumeric);
p.addRequired('Spin', @isnumeric);
p.addRequired('J_min', @isnumeric);
p.parse(B, L, N, Spin, J_min);
args = p.Results;

[psi_lm phi_l] = s2let_transform_tiling_mex(B, L, N, Spin, J_min);

end

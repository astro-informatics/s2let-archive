function [psi_lm phi_l] = s2let_transform_tiling(B, L, N, Spin, J_min, varargin)

% s2let_transform_tiling - Compute tiling in harmonic space.
% -- Directional wavelets on the sphere.
%
% Default usage :
%
%   [psi_lm phi_l] = s2let_transform_tiling(B, L, N, Spin, J_min, <options>)
%
% psi_lm is an array containing wavelet tiling.
% phi_l contains the scaling function.
% B is the wavelet parameter,
% L is the angular band-limit,
% N is the azimuthal band-limit,
% Spin is the spin number,
% J_min the first wavelet to be used.
%
% Options consist of parameter type and value pairs.
% Valid options include:
%
%  'SpinLowered' = { true  [Apply normalisation factors for spin-lowered
%                           wavelets and scaling function.],
%                    false [Apply the usual normalisation factors such
%                           that the wavelets fulfil the admissibility
%                           condition (default)]}
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2014 Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

p = inputParser;
p.addRequired('B', @isnumeric);
p.addRequired('L', @isnumeric);
p.addRequired('N', @isnumeric);
p.addRequired('Spin', @isnumeric);
p.addRequired('J_min', @isnumeric);
p.addParamValue('SpinLowered', false, @islogical);
p.parse(B, L, N, Spin, J_min, varargin{:});
args = p.Results;

[psi_lm phi_l] = s2let_transform_tiling_mex(B, L, N, Spin, J_min, args.SpinLowered);

end

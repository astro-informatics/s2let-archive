function f = s2let_compute_wav(j, L, alpha, beta, gamma, varargin)

% s2let_compute_wav 
% Compute the j-th wavelet at rho=(alpha, beta, gamma) on SO3
% and reconstruct it on the sphere
%
% Default usage :
%
%   f = s2let_compute_wav(alpha, beta, gamma, j, <options>)
%
% j is the order of the wavelet under consideration (depends on B and J_min)
% rho=(alpha, beta, gamma) is the position on SO3 where to evaluate the wavelet
% L if Harmonic band-limit for the reconstruction on the sphere
% f is the output reconstructed wavelet on the sphere, at resolution L
% 
%
% Option :
%  'B'               = { Dilation factor; B > 1 (default=2) }
%  'N'               = { Directional band-limit; N > 1 (default=L) }
%  'Spin'            = { Spin; Spin >= 0 (default=0) }
%
% S2LET package to perform Wavelet transform on the Sphere.
% Copyright (C) 2012-2014  Boris Leistedt, Martin Büttner & Jason McEwen
% See LICENSE.txt for license details


end
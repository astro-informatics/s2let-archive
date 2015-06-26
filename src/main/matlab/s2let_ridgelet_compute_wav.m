function [ridgelet_wav, ridgelet_scal] = ...
   s2let_ridgelet_compute_wav(L, varargin)
% s2let_ridgelet_compute_wav - Compute ridgelets
%
% Compute ridgelet wavelet and scaling functions.
%
% Default usage :
%
%   [ridgelet_wav, ridgelet_scal] = s2let_compute_wav(L, <options>)
%
% where L is harmonic band-limit for the reconstruction on the sphere,
% [ridgelet_wav, ridgelet_scal] are the reconstruct ridgelet wavelets and
% scaling functions on the sphere for each j.
%
% Options consist of parameter type and value pairs.
% Valid options include:
%
%  'B'               = { Dilation factor; B > 1 (default = 2) }
%  'Spin'            = { Spin number; Spin >= 0 (default = 0) }
%  'J_min'           = { Minimum wavelet scale to consider;
%                        0 <= J_min < log_B(L) (default=0) }
%  'Sampling'        = { 'MW'           [McEwen & Wiaux sampling (default)],
%                        'MWSS'         [McEwen & Wiaux symmetric sampling] }
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%
% S2LET package to perform Wavelet transform on the Sphere.
% Copyright (C) 2012-2014  Boris Leistedt, Martin Büttner & Jason McEwen
% See LICENSE.txt for license details

% Parse arguments.
p = inputParser;
p.addRequired('L', @isnumeric);
p.addParamValue('B', 2, @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('J_min', 0, @isnumeric);
p.addParamValue('Sampling', 'MW', @ischar);
p.addParamValue('Reality', false, @islogical);
p.parse(L, varargin{:});
args = p.Results;

N = 1;

% Compute ring in harmonic space.
ring_lm = zeros(args.L^2,1);
for el = 0:args.L-1
        
   elon2 = el./2.0;
   logp = gammaln(2*elon2+1) - 2*elon2 * log(2) - 2 * gammaln(elon2+1);
   p0 = real((-1).^elon2) .* exp(logp);
    
   ind = ssht_elm2ind(el, 0);
   ring_lm(ind) = 2*pi * sqrt((2*el+1)/(4*pi)) * p0;
       
end

% Compute ring on the sphere.
ring = ssht_inverse(ring_lm, L, ...
   'Reality', true, ...
   'Method', args.Sampling, ...
   'Spin', args.Spin);

% Compute ridgelets.
[ridgelet_wav, ridgelet_scal] = ...
   s2let_transform_analysis_mw(ring, 'B', args.B, 'J_min', args.J_min, ...
   'N', N, 'Upsample', true, 'Spin', args.Spin, ...
   'Reality', args.Reality, 'Sampling', args.Sampling);

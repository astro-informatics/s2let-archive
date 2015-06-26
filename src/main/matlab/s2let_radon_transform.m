function [f_radon_lm] = s2let_radon_transform(f_lm, varargin)

% s2let_radon_transform
% Compute radon transform in harmonic space.
%
% Default usage:
%
%   f_radon_lm = s2let_radon_transform(f_lm, <options>)
%
% where f_lm is the vector of L^2 harmonic coefficients and f_random_lm
% is the vector of harmonic coefficients of the spherical radon transform
% of f.
%
% Option :
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed from input) }

% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2015 Jason McEwen
% See LICENSE.txt for license details

L_guess = sqrt(length(f_lm));

p = inputParser;
p.addRequired('f_lm', @isnumeric);
p.addParamValue('L', L_guess, @isnumeric);
p.addParamValue('Reality', false, @islogical);
p.parse(f_lm, varargin{:});
args = p.Results;

% Compute radon transform
ring_lm = zeros(args.L^2,1);
f_radon_lm = zeros(args.L^2,1);
for el = 0:args.L-1
        
   elon2 = el./2.0;
   logp = gammaln(2*elon2+1) - 2*elon2 * log(2) - 2 * gammaln(elon2+1);
   p0 = real((-1).^elon2) .* exp(logp);
    
   ind = ssht_elm2ind(el, 0);
   ring_lm(ind) = 2*pi * sqrt((2*el+1)/(4*pi)) * p0;
       
   if args.Reality
      m_min = 0;
   else
      m_min = -el;
   end

   for m = m_min:el
      ind_lm = ssht_elm2ind(el, m);
      f_radon_lm(ind_lm) = sqrt(4 * pi ./ (2*el+1)) ...
         * f_lm(ind_lm) * ring_lm(ind);
   end
       
end


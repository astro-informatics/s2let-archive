function [f_lm] = s2let_radon_inverse(f_radon_lm, varargin)

% s2let_radon_inverse
% Compute inverse radon transform in harmonic space.
%
% Default usage:
%
%   f_lm = s2let_radon_analysis(f_radon_lm, <options>)
%
% where f_lm is the vector of L^2 harmonic coefficients and f_random_lm
% is the vector of harmonic coefficients of the spherical radon transform
% of f.
%
% Note that the inverse random transform recovers the antipodal part
% of the original signal f only.
%
% Option :
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed from input) }

% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2015 Jason McEwen
% See LICENSE.txt for license details

L_guess = sqrt(length(f_radon_lm));

p = inputParser;
p.addRequired('f_lm', @isnumeric);
p.addParamValue('L', L_guess, @isnumeric);
p.addParamValue('Reality', false, @islogical);
p.parse(f_radon_lm, varargin{:});
args = p.Results;

ring_lm = zeros(args.L^2,1);
f_lm = zeros(args.L^2,1);
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
      
      if mod(el, 2) == 1
         f_lm(ind_lm) = 0.0;
      else         
         f_lm(ind_lm) = f_radon_lm(ind_lm) ...
            ./ sqrt(4 * pi ./ (2*el+1)) ./ ring_lm(ind);
      end
   
   end
   
end

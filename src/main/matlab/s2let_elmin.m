function elmin = s2let_elmin(B, j)

% s2let_elmin
% Return the lowest harmonic index el supported by the given
% wavelet scale.
%
% Default usage:
%		elmin = s2let_elmin(B, j)
%
%   B is the wavelet parameter,
%   j the scale of interest.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

elmin = ceil(B^j);

end

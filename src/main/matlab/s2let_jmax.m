function J = s2let_jmax(L, B)

% s2let_jmax
% Return the maximum scale for a wavelet transform
%
% Default usage:
%		J_max = s2let_jmax(L, B)
$	
%	L is the band-limit for the transform,
%	B is the wavelet parameter.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

J = ceil(log(L) ./ log(B));

end

function f = s2let_axisym_synthesis(f_wav, f_scal, B, L, J_min, varargin)

% s2let_axisym_synthesis 
% Compute axisymmetric wavelet transform, output in pixel space.
%
% Default usage :
%
%   f = s2let_axisym_synthesis(f_wav, f_scal, B, L, J_min, <options>)
%
% f_wav contains the input wavelet contributions -- MW sampling,
% f_scal contains the input scaling contributions -- MW sampling,
% f is the output field -- MW sampling,
% B is the wavelet parameter,
% L is the angular band-limit,
% J_min the first wavelet to be used.
%
% Option :
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

p = inputParser;
p.addRequired('f_wav'); 
p.addRequired('f_scal', @isnumeric); 
p.addRequired('B', @isnumeric);   
p.addRequired('L', @isnumeric);   
p.addRequired('J_min', @isnumeric); 
p.addParamValue('Reality', false, @islogical);
p.parse(f_wav, f_scal, B, L, J_min, varargin{:});
args = p.Results;

f_scal_vec = s2let_mw_arr2vec(f_scal);

J = s2let_jmax(L, B);
f_wav_vec = zeros(J+1, L*(2*L-1));
for j = 0:J
	temp = f_wav{j+1};
  	f_wav_vec(j+1,:) = s2let_mw_arr2vec(temp);
end

f_vec = s2let_axisym_synthesis_mex(f_wav_vec, f_scal_vec, B, L, J_min, args.Reality);
 
f = s2let_mw_vec2arr(f_vec);

end
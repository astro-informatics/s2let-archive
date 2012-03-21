function [f_wav, f_scal] = s2let_axisym_analysis(f, B, L, J_min, varargin)

% s2let_axisym_analysis 
% Compute axisymmetric wavelet transform, output in pixel space.
%
% Default usage :
%
%   [f_wav, f_scal] = s2let_axisym_analysis(f, B, L, J_min, <options>)
%
% f is the input field -- MW sampling,
% f_wav contains the output wavelet contributions,
% f_scal contains the output scaling contributions,
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
p.addRequired('f', @isnumeric); 
p.addRequired('B', @isnumeric);   
p.addRequired('L', @isnumeric);   
p.addRequired('J_min', @isnumeric); 
p.addParamValue('Reality', false, @islogical);
p.parse(f, B, L, J_min, varargin{:});
args = p.Results;

f_vec = s2let_mw_arr2vec(f);

[f_wav_vec, f_scal_vec] = s2let_axisym_analysis_mex(f_vec, B, L, J_min, args.Reality);
 
f_scal = s2let_mw_vec2arr(f_scal_vec);

J = s2let_jmax(L, B);
f_wav = cell(J+1, 1);
for j = 0:J
  f_wav{j+1} = s2let_mw_vec2arr(f_wav_vec(j+1,:));
end

end
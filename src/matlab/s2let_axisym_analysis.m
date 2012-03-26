function [f_wav, f_scal] = s2let_axisym_analysis(f, varargin)

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
%
% Option :
%  'B'               = { Dilation factor; B > 1 (default=2) }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed) }
%  'J_min'           = { Minimum needlet scale to consider;
%                        0 <= J_min < log_B(L) (default=0) }
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(f);
Lguessed = min([sz(1) sz(2)]);

p = inputParser;
p.addRequired('f', @isnumeric); 
p.addParamValue('B', 2, @isnumeric);   
p.addParamValue('L', Lguessed, @isnumeric);   
p.addParamValue('J_min', 0, @isnumeric); 
p.addParamValue('Reality', false, @islogical);
p.parse(f, varargin{:});
args = p.Results;

f_vec = s2let_mw_arr2vec(f);

[f_wav_vec, f_scal_vec] = s2let_axisym_analysis_mex(f_vec, args.B, args.L, args.J_min, args.Reality);
 
f_scal = s2let_mw_vec2arr(f_scal_vec);

J = s2let_jmax(args.L, args.B);
f_wav = cell(J+1, 1);
for j = 0:J
  f_wav{j+1} = s2let_mw_vec2arr(f_wav_vec(j+1,:));
end

end
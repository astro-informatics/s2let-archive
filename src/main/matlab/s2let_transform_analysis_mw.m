function [f_wav, f_scal] = s2let_transform_analysis_mw(f, varargin)

% s2let_transform_analysis_mw
% Compute spin directional wavelet transform, output in pixel space.
%
% Default usage :
%
%   [f_wav, f_scal] = s2let_transform_analysis_mw(f, <options>)
%
% f is the input field -- MW sampling,
% f_wav contains the output wavelet contributions,
% f_scal contains the output scaling contributions,
%
% Option :
%  'B'               = { Dilation factor; B > 1 (default=2) }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed from input) }
%  'N'               = { Azimuthal/directional band-limit; N > 1 (default=L) }
%  'Spin'               = { Spin; (default=0) }
%  'J_min'           = { Minimum needlet scale to consider;
%                        0 <= J_min < log_B(L) (default=0) }
%  'Downsample'      = { true        [multiresolution algorithm (default)],
%                        false       [full resolution wavelets] }
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%  'SpinLowered'     = { true  [Apply normalisation factors for spin-lowered
%                               wavelets and scaling function.],
%                        false [Apply the usual normalisation factors such
%                               that the wavelets fulfil the admissibility
%                               condition (default)]}
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
p.addParamValue('N', Lguessed, @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Downsample', true, @islogical);
p.addParamValue('Reality', false, @islogical);
p.addParamValue('SpinLowered', false, @islogical)
p.parse(f, varargin{:});
args = p.Results;

f_vec = s2let_mw_arr2vec(f);

[f_wav_vec, f_scal_vec] = s2let_transform_analysis_mw_mex(f_vec, args.B, args.L, args.J_min, args.N, args.Spin, args.Reality, args.Downsample, args.SpinLowered);

f_scal = s2let_mw_vec2arr(f_scal_vec);

J = s2let_jmax(args.L, args.B);
f_wav = cell(J+1-args.J_min, 2*args.N-1);
offset = 0;
for j = args.J_min:J
  for en = 1:2*args.N-1
    if args.Downsample == true
      band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
    else
      band_limit = args.L;
    end
    temp = zeros(band_limit, 2*band_limit-1);
    for t = 0:band_limit-1
        for p = 0:2*band_limit-2
            ind = offset + t * ( 2 * band_limit - 1) + p + 1;
            temp(t+1,p+1) = f_wav_vec(1,ind);
        end
    end
    f_wav{j+1-args.J_min, en} = temp;
    offset = offset + band_limit * (2*band_limit-1);
  end
end

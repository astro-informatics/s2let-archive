
% s2let_demo_ridgelet_plot
% Compute and plot the ridgelet wavelet and scaling functions.
%
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

B = 2;
N = 1;
J_min = 3;
J = s2let_jmax(L, B);

j_plot = J_min + 1;
[ridgelet_wav, ridgelet_scal] = s2let_ridgelet_compute_wav(L, ...
   'B', B, 'J_min', J_min, ...
   'Spin', 0, 'Reality', true, 'Sampling', sampling_method);

figure
ssht_plot_sphere(ridgelet_wav{j_plot-J_min+1}, L, ...
   'Lighting', sphere_plot_lighting, ...
   'Type', sphere_plot_type, 'ParametricScale', sphere_plot_scale, ...
   'Method', sampling_method);

figure
ssht_plot_sphere(ridgelet_scal, L, ...
   'Lighting', sphere_plot_lighting, ...
   'Type', sphere_plot_type, 'ParametricScale', sphere_plot_scale, ...
   'Method', sampling_method);

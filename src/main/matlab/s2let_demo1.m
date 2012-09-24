% s2let_demo1
% Analyse Earth tomography data as a real MW map.
% Compute the wavelet maps and plot them.
% Plot 1 : multiresolution wavelet scales
% Plot 2 : full resolution wavelet scales
%
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

load('EGM2008_Topography_flms_L0128');
f = ssht_inverse(flm, L, 'Reality', true);

%inputfile = 'data/earth_tomo_mw_128.fits';
%[f, L] = s2let_read_mw_real_map(inputfile);

B = 3;
J_min = 2;
J = s2let_jmax(L, B);

% Perform decomposition
[f_wav, f_scal] = s2let_axisym_analysis(f, 'B', B, 'J_min', J_min, 'Reality', true);

zoomfactor = 1.2;
ns = ceil(sqrt(2+J-J_min+1)) ;
ny = ns - 1 + rem(2+J-J_min + 1, ns) ;
nx = ns;
% MULTIRESOLUTION PLOT
figure('Position',[100 100 1300 1000])
subplot(nx, ny, 1);
s2let_plot_mw_mollweide(f);
%title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
subplot(nx, ny, 2);
s2let_plot_mw_mollweide(f_scal);
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
%title('Scaling fct')
for j = J_min:J
   subplot(nx, ny, j-J_min+3);
   s2let_plot_mw_mollweide(f_wav{j-J_min+1});
   campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
   %title(['Wavelet scale : ',int2str(j)-J_min+1])
end 

% Perform decomposition
[f_wav, f_scal] = s2let_axisym_analysis(f, 'B', B, 'J_min', J_min, 'Reality', true, 'downsample', false);
% FULL RESOLUTION PLOT
figure('Position',[100 100 1300 1000])
subplot(nx, ny, 1);
s2let_plot_mw_mollweide(f);
%title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
subplot(nx, ny, 2);
s2let_plot_mw_mollweide(f_scal);
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
%title('Scaling fct')
for j = J_min:J
   subplot(nx, ny, j-J_min+3);
   s2let_plot_mw_mollweide(f_wav{j-J_min+1});
   campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
   %title(['Wavelet scale : ',int2str(j)-J_min+1])
end

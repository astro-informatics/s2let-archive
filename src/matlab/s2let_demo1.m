% s2let_demo1
% Analyse Earth tomography data
% Plot 1 : multiresolution wavelet scales
% Plot 2 : full resolution wavelet scales

load('EGM2008_Topography_flms_L0128');
f = ssht_inverse(flm, L, 'Reality', true);

B = 2;
J_min = 0;
J = s2let_jmax(L, B);

% Perform decomposition
[f_wav, f_scal] = s2let_axisym_analysis(f, 'B', B, 'J_min', J_min, 'Reality', true);

% MULTIRESOLUTION PLOT
figure('Position',[100 100 1400 1000])
subplot((J-J_min+1)/2-1, (J-J_min+1)/2-1, 1);
ssht_plot_mollweide(f, L);
title('Initial data')
zoom(1.5)
subplot((J-J_min+1)/2-1, (J-J_min+1)/2-1, 2);
band_limit = min([ ceil(B^(J_min+1)) L ]);
ssht_plot_mollweide(f_scal, band_limit);
zoom(1.5)
title('Scaling fct')
for j = J_min:J
   subplot((J-J_min+1)/2-1, (J-J_min+1)/2-1, j-J_min+2);
   band_limit = min([ ceil(B^(j+1)) L ]);
   ssht_plot_mollweide(f_wav{j-J_min+1}, band_limit);
   zoom(1.5)
   title(['Wavelet cale : ',int2str(j)])
end 

% Perform decomposition
[f_wav, f_scal] = s2let_axisym_analysis(f, 'B', B, 'J_min', J_min, 'Reality', true, 'downsample', false);

% FULL RESOLUTION PLOT
figure('Position',[100 100 1400 1000])
subplot((J-J_min+1)/2-1, (J-J_min+1)/2-1, 1);
ssht_plot_mollweide(f, L);
title('Initial data')
zoom(1.5)
subplot((J-J_min+1)/2-1, (J-J_min+1)/2-1, 2);
ssht_plot_mollweide(f_scal, L);
zoom(1.5)
title('Scaling fct')
for j = J_min:J
   subplot((J-J_min+1)/2-1, (J-J_min+1)/2-1, j-J_min+2);
   ssht_plot_mollweide(f_wav{j-J_min+1}, L);
   zoom(1.5)
   title(['Wavelet scale : ',int2str(j)])
end

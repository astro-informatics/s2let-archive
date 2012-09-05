% s2let_demo1
% Analyse Earth tomography data
% Plot 1 : multiresolution wavelet scales
% Plot 2 : full resolution wavelet scales

load('EGM2008_Topography_flms_L0128');
f = ssht_inverse(flm, L, 'Reality', true);

% load('wmap_lcdm_pl_model_wmap7baoh0.mat');
% L = 256;
% cmb_lm = zeros(L^2,1);
% for el = 2:L-1      
%    cl(el-1) = cl(el-1) * 2*pi / (el*(el+1));
%    for m = -el:el
%      ind = ssht_elm2ind(el, m);
%      if m == 0
%        cmb_lm(ind) = sqrt(cl(el-1)) .* randn;
%      else
%         cmb_lm(ind) = ...
%            sqrt(cl(el-1)./2) .* randn ...
%            + sqrt(-1) * sqrt(cl(el-1)./2) .* randn;
%      end
%    end
% end
% 
% % Compute real space version of CMB.
% f = ssht_inverse(cmb_lm, L, 'Reality', true);

B = 3;
J_min = 2;
J = s2let_jmax(L, B)

% Perform decomposition
[f_wav, f_scal] = s2let_axisym_analysis(f, 'B', B, 'J_min', J_min, 'Reality', true);

zoomfactor = 1.2;
ns = ceil(sqrt(2+J-J_min+1)) ;
ny = ns - 1 + rem(2+J-J_min + 1, ns) ;
nx = ns;

% MULTIRESOLUTION PLOT
figure('Position',[100 100 1300 1000])
subplot(nx, ny, 1);
ssht_plot_mollweide(f, L);
title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
subplot(nx, ny, 2);
band_limit = min([ ceil(B^(J_min)) L ]);
ssht_plot_mollweide(f_scal, band_limit);
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
title('Scaling fct')
for j = J_min:J
   subplot(nx, ny, j-J_min+3);
   band_limit = min([ ceil(B^(j+1)) L ]);
   ssht_plot_mollweide(f_wav{j-J_min+1}, band_limit);
   campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
   title(['Wavelet scale : ',int2str(j)-J_min+1])
end 

% Perform decomposition
[f_wav, f_scal] = s2let_axisym_analysis(f, 'B', B, 'J_min', J_min, 'Reality', true, 'downsample', false);

% FULL RESOLUTION PLOT
figure('Position',[100 100 1300 1000])
subplot(nx, ny, 1);
ssht_plot_mollweide(f, L);
title('Initial data')
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
subplot(nx, ny, 2);
ssht_plot_mollweide(f_scal, L);
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
title('Scaling fct')
for j = J_min:J
   subplot(nx, ny, j-J_min+3);
   ssht_plot_mollweide(f_wav{j-J_min+1}, L);
   campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
   title(['Wavelet scale : ',int2str(j)-J_min+1])
end

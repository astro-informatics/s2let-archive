% s2let_demo3
% Analyse some CMB simulation from a valid HEALPIX map

L = 128;
B = 3;
J_min = 1;

[f, nside] = read_healpix_map('data/some_cmb_simu.fits');

% Perform decomposition
[f_wav, f_scal] = s2let_hpx_axisym_analysis(f,'B',B,'L',L,'J_min',J_min);

% Perform reconstruction
f_rec = s2let_hpx_axisym_synthesis(f_wav, f_scal,'B',B,'L',L,'J_min',J_min);


% Plot

J = s2let_jmax(L, B);
zoomfactor = 1.2;
ns = ceil(sqrt(2+J-J_min+1)) ;
ny = ns - 1 + rem(2+J-J_min + 1, ns) ;
nx = ns;
figure('Position',[100 100 1300 1000])

subplot(nx, ny, 1);
s2let_plot_hpx_mollweide(f);
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
title('Initial data')

subplot(nx, ny, 2);
s2let_plot_hpx_mollweide(f_rec);
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
title('Reconstruction')

subplot(nx, ny, 3);
s2let_plot_hpx_mollweide(f_scal);
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
title('Scaling fct')

for j = J_min:J
   subplot(nx, ny, j-J_min+4);
   s2let_plot_hpx_mollweide(f_wav{j-J_min+1});
   campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
   title(['Wavelet scale : ',int2str(j)-J_min+1])
end 

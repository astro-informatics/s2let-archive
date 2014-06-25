%% Preparation

position = [ 500 500 1000 800 ];
set(0, 'DefaultFigurePosition', position);
plot_root = '../../../figs/';
plot_size = [0 0 23 13];

% Set working directory to script's location
% (assumed to be src/main/matlab)

cd(fileparts(mfilename('fullpath')))

%% Plot real map results

[f, L] = s2let_mw_read_real_map('../../../data/real_signal_input.fits');
[f_noisy, ~] = s2let_mw_read_real_map('../../../data/real_signal_noisy.fits');
[f_denoised, ~] = s2let_mw_read_real_map('../../../data/real_signal_denoised.fits');

% Precompute Wigner small-d functions
sqrt_tbl = sqrt([0:2*(L-1)+1])';
signs = ones(L+1,1);
signs(2:2:end) = -1; 

d = zeros(L, 2*L-1, 2*L-1);
d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), L, 0, 0, ...
  'SqrtTable', sqrt_tbl, 'SignTable', signs);
for el = 1:L-1
  d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), L, el, 0, ...
     'SqrtTable', sqrt_tbl, 'SignTable', signs);
end

% Rotate maps by 180°
f = ssht_inverse(ssht_rotate_flm(ssht_forward(f, L, 'Reality', true), d, pi, 0), L, 'Reality', true);
f_noisy = ssht_inverse(ssht_rotate_flm(ssht_forward(f_noisy, L, 'Reality', true), d, pi, 0), L, 'Reality', true);
f_denoised = ssht_inverse(ssht_rotate_flm(ssht_forward(f_denoised, L, 'Reality', true), d, pi, 0), L, 'Reality', true);

f_error = log(abs(f_denoised./f - 1));

figure
ssht_plot_mollweide(f, L, 'Mode', 1)
colorbar;
title('Real signal - Input map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_real_signal_input.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_noisy, L, 'Mode', 1)
colorbar;
title('Real signal - Input map with added noise');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_real_signal_noisy.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_denoised, L, 'Mode', 1)
colorbar;
title('Real signal - Denoised map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_real_signal_denoised.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_error, L, 'Mode', 1)
colorbar;
title('Real signal - Logarithmic relative error');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_real_signal_error.png'];
print('-r200', '-dpng', fname)

%% Plot spin map results

[f, L] = s2let_mw_read_real_map('../../../data/spin_signal_real_input.fits');
[f_noisy, ~] = s2let_mw_read_real_map('../../../data/spin_signal_real_noisy.fits');
[f_denoised, ~] = s2let_mw_read_real_map('../../../data/spin_signal_real_denoised.fits');

% Precompute Wigner small-d functions
sqrt_tbl = sqrt([0:2*(L-1)+1])';
signs = ones(L+1,1);
signs(2:2:end) = -1; 

d = zeros(L, 2*L-1, 2*L-1);
d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), L, 0, 0, ...
  'SqrtTable', sqrt_tbl, 'SignTable', signs);
for el = 1:L-1
  d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), L, el, 0, ...
     'SqrtTable', sqrt_tbl, 'SignTable', signs);
end

% Rotate maps by 180°
f = ssht_inverse(ssht_rotate_flm(ssht_forward(f, L, 'Reality', true), d, pi, 0), L, 'Reality', true);
f_noisy = ssht_inverse(ssht_rotate_flm(ssht_forward(f_noisy, L, 'Reality', true), d, pi, 0), L, 'Reality', true);
f_denoised = ssht_inverse(ssht_rotate_flm(ssht_forward(f_denoised, L, 'Reality', true), d, pi, 0), L, 'Reality', true);

f_error = log(abs(f_denoised./f - 1));

figure
ssht_plot_mollweide(f, L, 'Mode', 1)
colorbar;
title('Spin signal (Q) - Input map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_Q_input.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_noisy, L, 'Mode', 1)
colorbar;
title('Spin signal (Q) - Input map with added noise');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_Q_noisy.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_denoised, L, 'Mode', 1)
colorbar;
title('Spin signal (Q) - Denoised map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_Q_denoised.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_error, L, 'Mode', 1)
colorbar;
title('Spin signal (Q) - Logarithmic relative error');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_Q_error.png'];
print('-r200', '-dpng', fname)

[f, ~] = s2let_mw_read_real_map('../../../data/spin_signal_imag_input.fits');
[f_noisy, ~] = s2let_mw_read_real_map('../../../data/spin_signal_imag_noisy.fits');
[f_denoised, ~] = s2let_mw_read_real_map('../../../data/spin_signal_imag_denoised.fits');

% Rotate maps by 180°
f = ssht_inverse(ssht_rotate_flm(ssht_forward(f, L, 'Reality', true), d, pi, 0), L, 'Reality', true);
f_noisy = ssht_inverse(ssht_rotate_flm(ssht_forward(f_noisy, L, 'Reality', true), d, pi, 0), L, 'Reality', true);
f_denoised = ssht_inverse(ssht_rotate_flm(ssht_forward(f_denoised, L, 'Reality', true), d, pi, 0), L, 'Reality', true);

f_error = log(abs(f_denoised./f - 1));

figure
ssht_plot_mollweide(f, L, 'Mode', 1)
colorbar;
title('Spin signal (U) - Input map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_U_input.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_noisy, L, 'Mode', 1)
colorbar;
title('Spin signal (U) - Input map with added noise');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_U_noisy.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_denoised, L, 'Mode', 1)
colorbar;
title('Spin signal (U) - Denoised map');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_U_denoised.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_error, L, 'Mode', 1)
colorbar;
title('Spin signal (U) - Logarithmic relative error');
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_U_error.png'];
print('-r200', '-dpng', fname)

%% Clean up

close all

clear position plot_root plot_size
clear f f_noisy f_denoised f_error
clear L
clear fname

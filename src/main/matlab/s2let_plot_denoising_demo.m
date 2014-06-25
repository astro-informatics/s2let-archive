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
f_error = abs(f_denoised./f - 1);

figure
ssht_plot_mollweide(f, L, 'Mode', 1)

set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_real_signal_input.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_noisy, L, 'Mode', 1)
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_real_signal_noisy.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_denoised, L, 'Mode', 1)
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_real_signal_denoised.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_error, L, 'Mode', 1)
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_real_signal_error.png'];
print('-r200', '-dpng', fname)

%% Plot spin map results

[f, L] = s2let_mw_read_real_map('../../../data/spin_signal_real_input.fits');
[f_noisy, ~] = s2let_mw_read_real_map('../../../data/spin_signal_real_noisy.fits');
[f_denoised, ~] = s2let_mw_read_real_map('../../../data/spin_signal_real_denoised.fits');
f_error = abs(f_denoised./f - 1);

figure
ssht_plot_mollweide(f, L, 'Mode', 1)

set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_Q_input.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_noisy, L, 'Mode', 1)
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_Q_noisy.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_denoised, L, 'Mode', 1)
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_Q_denoised.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_error, L, 'Mode', 1)
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_Q_error.png'];
print('-r200', '-dpng', fname)

[f, L] = s2let_mw_read_real_map('../../../data/spin_signal_imag_input.fits');
[f_noisy, ~] = s2let_mw_read_real_map('../../../data/spin_signal_imag_noisy.fits');
[f_denoised, ~] = s2let_mw_read_real_map('../../../data/spin_signal_imag_denoised.fits');
f_error = abs(f_denoised./f - 1);

figure
ssht_plot_mollweide(f, L, 'Mode', 1)

set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_U_input.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_noisy, L, 'Mode', 1)
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_U_noisy.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_denoised, L, 'Mode', 1)
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', plot_size);
colormap(jet);
fname = [plot_root, 'denoising_spin_signal_U_denoised.png'];
print('-r200', '-dpng', fname)

figure
ssht_plot_mollweide(f_error, L, 'Mode', 1)
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

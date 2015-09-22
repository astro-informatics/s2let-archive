make mw_bin
bin/s2let_denoising_demo
bin/s2let_spin_denoising_demo

matlab -nodesktop -nosplash -r "s2let_plot_denoising_demo;exit"

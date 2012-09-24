function s2let_plot_mw_mollweide(f)
% s2let_plot_hpx_mollweide - Plot HEALPIX map using Mollweide projection

sz = size(f);
L = min(sz);
ssht_plot_mollweide(f, L);
campos([0 0 -1]); camup([0 1 0]);

end
function s2let_plot_mw_mollweide(f)

% s2let_plot_mw_mollweide 
% Plot a real MW map using Mollweide projection.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(f);
L = min(sz);
ssht_plot_mollweide(f, L);
campos([0 0 -1]); camup([0 1 0]);

end
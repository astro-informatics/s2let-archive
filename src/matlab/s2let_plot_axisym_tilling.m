function s2let_plot_axisym_tilling(B, L, J_min)

% plot_axisym_tilling - Plot tilling in harmonic space.
% -- Axisymmetric wavelets on the sphere.
%
% Default usage :
%
%   plot_axisym_tilling(B, L, J_min)
%
% B is the wavelet parameter,
% L is the angular band-limit,
% J_min the first wavelet to be used.
%
% S2LET package to perform Wavelet transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

[kappa kappa0] = s2let_axisym_tilling(B, L, J_min);

J = s2let_jmax(L, B);

figure;
plot(0:L-1, kappa0, 'k', 'LineWidth', 2);
h = text(1, 1.1, 'kappa0', 'Color', [0 0 0]);
hold on;
for j = J_min:J  
  colour = rand(1,3)*0.9;
  plot(0:L-1, kappa(j+1,:), 'LineWidth', 2, 'Color', colour);
  h = text(B.^j, 1.05, strcat('j=',num2str(j)), 'Color', colour);  
end
title('Harmonic tiling');
xlabel('el');
axis([0 L -0.05 1.15]);
set(gca,'XTick',B.^[0:J]);

end
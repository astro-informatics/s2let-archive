% s2let_demo2
% plot wavelet kernels

B = 4;
zoomfactor = 1.2;
J_min = 1;
L = 256;
J = s2let_jmax(L, B);
type = 'colour';
lighting = true;

ns = ceil(sqrt(2+J-J_min+1)) ;
ny = ns - 1 + rem(2+J-J_min + 1, ns) ;
nx = ns;

[kappa kappa0] = s2let_axisym_tilling(B, L, J_min);

s2let_plot_axisym_tilling(B, L, J_min);

figure('Position',[100 100 1100 700]) 

subplot(nx, ny, 1);
flm = zeros(L^2,1);
for l = 0:L-1
    flm(l^2+l+1,1) = kappa0(l+1);
end     
f = ssht_inverse(flm, L, 'Reality', true);
ssht_plot_sphere(f, L, 'Type', type, 'Lighting', lighting);
zoom(zoomfactor)
%title('Scaling fct')

for j = J_min:J
   subplot(nx, ny, j-J_min+2);
   flm = zeros(L^2,1);
    for l = 0:L-1
        flm(l^2+l+1,1) = kappa(j+1,l+1);
    end  
   f = ssht_inverse(flm, L, 'Reality', true);
   ssht_plot_sphere(f, L, 'Type', type, 'Lighting', lighting);
   %title(['Wavelet scale : ',int2str(j)-J_min+1])
   zoom(zoomfactor)
end
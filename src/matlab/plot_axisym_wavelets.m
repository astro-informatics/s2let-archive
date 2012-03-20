function plot_laguerre_kernels(N, R, npoints)

h = R/(npoints);
nodes = (h:h:R);
tau = flag_get_tau(N, R);
sampling = slag_sampling(N, R);


            
nmax = N;%10;
colors = colorrange(nmax);

figure('Position',[1 1 500 1000])

subplot(2,1,1)
hold on
for n = 1:nmax
    fn = zeros(1,N);
    fn(n) = 1.0;
    [f, nodes] = slag_synthesis(fn, N, 'Nodes', nodes);
    p = plot(nodes, f.*nodes/sqrt(tau));
    set(p,'Color',colors(n,:),'LineWidth',1.1);
end
plot(sampling, 0*sampling,'d','LineWidth',1.1,'MarkerSize',10,'MarkerEdgeColor','k',...
                'MarkerFaceColor',[.49 1 .63])
axis([0 R -1 1])
xlabel('Radius')
title('Decaying Laguerre polynomials (weighted with exp-r)')

subplot(2,1,2)
hold on
for n = 1:nmax
    fn = zeros(1,N);
    fn(n) = 1.0;
    [f, nodes] = slag_synthesis(fn, N, 'Nodes', nodes);
    p = plot(nodes, f/sqrt(tau));
    set(p,'Color',colors(n,:),'LineWidth',1.1);
end
plot(sampling, 0*sampling,'d','LineWidth',1.1,'MarkerSize',10,'MarkerEdgeColor','k',...
                'MarkerFaceColor',[.49 1 .63])
axis([0 R -2 2])
xlabel('Radius')
title('Spherical Laguerre Basis functions')

end
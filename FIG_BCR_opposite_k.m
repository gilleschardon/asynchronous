clear all

Nm = 5;
Nm2 =  Nm^2;
K = [linspace(0.1, 0.5, 100) (1:50)/2];
Nsnaptot = 100;

sigma2 = 0.5;

xx = linspace(-0.25, 0.25, Nm);

[X, Y] = meshgrid(xx, xx);

Z1 = -2;
Z2 = 2;

XYZm1 = [X(:) Y(:) Z1 * ones(size(X(:)))];
XYZm2 = [X(:) Y(:) Z2 * ones(size(X(:)))];
%XYZm2 = [ Z2 * ones(size(X(:))) Y(:) X(:)];

XYZs = [0.1, 0.0, 0.2];

p = 1;


B3r = zeros(size(K));
B3s = zeros(size(K));
B3 = zeros(size(K));

Bps = zeros(size(K));
Bpr = zeros(size(K));
Bp = zeros(size(K));


for u = 1:length(K)
    k = K(u);
    [B3s(u), Bps(u), B3r(u), Bpr(u)] = BCRuncN({XYZm1, XYZm2}, [Nsnaptot, Nsnaptot], XYZs, k,  p, sigma2, @freefieldsource);
    [B3(u), Bp(u)] = BCRuncN({[XYZm1; XYZm2]}, [Nsnaptot+Nsnaptot], XYZs, k,  p, sigma2, @freefieldsource);


end
%%

figure
subplot(2, 2, 2)
plot(K, B3s./B3, '-', 'linewidth', 2)
hold on
%plot(K, B3r./B3s, '--', 'linewidth', 2)
plot(K, B3r./B3s, '--', 'linewidth', 2)


xlabel('k (m^{-1})')
ylabel('CRB ratio, position')
legend('$B_{\mathbf x}^s/ B_{\mathbf x}^c$', '$B_{\mathbf x}^r/ B_{\mathbf x}^s$' , 'fontsize', 15, 'interpreter', 'latex')



subplot(2, 2, 4)

plot(K, Bps./Bp, '-', 'linewidth', 2)
hold on
plot(K, Bpr./Bps, '--', 'linewidth', 2)

xlabel('k (m^{-1})')
ylabel('CRB ratio, power')
legend('$B_{\mathbf x}^s/ B_{\mathbf x}^c$', '$B_{\mathbf x}^r/ B_{\mathbf x}^s$' , 'fontsize', 15, 'interpreter', 'latex')

%,'$B_{\mathbf x}^r/ B_{\mathbf x}^s$','$B_{\mathbf x}^r/ B_{\mathbf x}^c$' , 'fontsize', 20, 'interpreter', 'latex')


subplot(2, 2, 1)
semilogy(K, B3, '-', 'linewidth', 2)
hold on
plot(K, B3s, '--', 'linewidth', 2)
plot(K, B3r, '-.', 'linewidth', 2)


xlabel('k (m^{-1})')
ylabel('CRB position (m^2)')
legend('$B_{\mathbf x}^c$', '$B_{\mathbf x}^s$', '$B_{\mathbf x}^r$' , 'fontsize', 15, 'interpreter', 'latex')



subplot(2, 2, 3)

plot(K, Bp, '-', 'linewidth', 2)
hold on
plot(K, Bps, '--', 'linewidth', 2)
plot(K, Bpr, '-.', 'linewidth', 2)

xlabel('k (m^{-1})')
ylabel('CRB power (Pa^2)')
legend('$B_{p}^c$', '$B_{p}^s$', '$B_{p}^r$' , 'fontsize', 15, 'interpreter', 'latex')

%,'$B_{\mathbf x}^r/ B_{\mathbf x}^s$','$B_{\mathbf x}^r/ B_{\mathbf x}^c$' , 'fontsize', 20, 'interpreter', 'latex')

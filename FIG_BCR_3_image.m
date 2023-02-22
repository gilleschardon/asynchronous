%clear all

Nm = 4;
Nm2 =  Nm^2;
k = 10;
Nsnaps = [100, 100, 100];

sigma2 = 0.2;

% Arrays coordinates

xx = linspace(-0.25, 0.25, Nm);

[X, Y] = meshgrid(xx, xx);

Z1 = -2;
Z2 = 2;

XYZm1 = [X(:) Y(:) Z1 * ones(size(X(:)))];

theta = 2*pi/3;
Rot = [cos(theta) , 0, sin(theta) ; 0, 1, 0 ; -sin(theta), 0, cos(theta)];
XYZm2 = [X(:) Y(:) Z1 * ones(size(X(:)))] * Rot;
XYZm3 = [X(:) Y(:) Z1 * ones(size(X(:)))] * Rot * Rot;


% initialization grid

resolg = 51;
xg = linspace(-1.8, 1.8, resolg);
yg = 0;
zg = linspace(-1.8, 1.8, resolg);

[XG, YG, ZG] = meshgrid(xg, yg, zg);
Xg = [XG(:) YG(:) ZG(:)];

p = 1;

% CRBs

B3r = zeros(size(Xg, 1), 1);
B3s = zeros(size(Xg, 1), 1);
B3 = zeros(size(Xg, 1), 1);

Bps = zeros(size(Xg, 1), 1);
Bpr = zeros(size(Xg, 1), 1);
Bp = zeros(size(Xg, 1), 1);

for u = 1:size(Xg, 1)

    [B3s(u), Bps(u), B3r(u), Bpr(u)] = BCRuncN({XYZm1, XYZm2, XYZm3}, Nsnaps, Xg(u, :), k,  p, sigma2, @freefieldsource);
    [B3(u), Bp(u)] = BCRuncN({[XYZm1; XYZm2; XYZm3]}, sum(Nsnaps), Xg(u, :), k,  p, sigma2, @freefieldsource);


end
%%

figure('Position', [100, 100, 1000, 600])


XLIM = [-2.1 2];
YLIM = [-2.1 2];
mask = (XG(:).^2 + ZG(:).^2) < 3.5;

subplot(2, 3, 1)

scatter(XYZm1(:, 1), XYZm1(:, 3), 'ko')
hold on
scatter(XYZm2(:, 1), XYZm2(:, 3), 'kx')
scatter(XYZm3(:, 1), XYZm3(:, 3), 'k+')

h = imagesc(xg, zg, reshape(B3, resolg, resolg)');

set(h, 'AlphaData', reshape(mask, resolg, resolg))
colormap(hot)
colorbar
xlabel('X (m)')
ylabel('Y (m)')
title('$B_{\mathbf x}^c$ (m$^2$)', 'interpreter', 'latex')
axis square

xlim(XLIM)
ylim(YLIM)

subplot(2, 3, 2)

scatter(XYZm1(:, 1), XYZm1(:, 3), 'ko')
hold on
scatter(XYZm2(:, 1), XYZm2(:, 3), 'kx')
scatter(XYZm3(:, 1), XYZm3(:, 3), 'k+')

h = imagesc(xg, zg, reshape(B3s./B3, resolg, resolg)');

set(h, 'AlphaData', reshape(mask, resolg, resolg))
colormap(hot)
colorbar
xlabel('X (m)')
ylabel('Y (m)')
title('$B_{\mathbf x}^s/ B_{\mathbf x}^c$', 'interpreter', 'latex')
caxis([1, max(B3s./B3.*mask)])

xlim(XLIM)
ylim(YLIM)
axis square

subplot(2, 3, 3)
scatter(XYZm1(:, 1), XYZm1(:, 3), 'ko')
hold on
scatter(XYZm2(:, 1), XYZm2(:, 3), 'kx')
scatter(XYZm3(:, 1), XYZm3(:, 3), 'k+')

h = imagesc(xg, zg, reshape(B3r./B3s, resolg, resolg)');
set(h, 'AlphaData', reshape(mask, resolg, resolg))

colormap(hot)
colorbar
xlabel('X (m)')
ylabel('Y (m)')
title('$B_{\mathbf x}^r/ B_{\mathbf x}^s$', 'interpreter', 'latex')
caxis([1, max(B3r./B3s.*mask)])
xlim(XLIM)
ylim(YLIM)
axis square


subplot(2, 3, 4)

scatter(XYZm1(:, 1), XYZm1(:, 3), 'ko')
hold on
scatter(XYZm2(:, 1), XYZm2(:, 3), 'kx')
scatter(XYZm3(:, 1), XYZm3(:, 3), 'k+')

h = imagesc(xg, zg, reshape(Bp, resolg, resolg)');

set(h, 'AlphaData', reshape(mask, resolg, resolg))
colormap(hot)
colorbar
xlabel('X (m)')
ylabel('Y (m)')
title('$B_{p}^c$ (Pa$^2$)', 'interpreter', 'latex')

xlim(XLIM)
ylim(YLIM)
axis square


subplot(2, 3, 5)
scatter(XYZm1(:, 1), XYZm1(:, 3), 'ko')
hold on
scatter(XYZm2(:, 1), XYZm2(:, 3), 'kx')
scatter(XYZm3(:, 1), XYZm3(:, 3), 'k+')

h = imagesc(xg, zg, reshape(Bps./Bp, resolg, resolg)');
set(h, 'AlphaData', reshape(mask, resolg, resolg))

colormap(hot)
colorbar
xlabel('X (m)')
ylabel('Y (m)')
title('$B_p^s/ B_p^c$', 'interpreter', 'latex')
caxis([1, max(Bps./Bp.*mask)])
xlim(XLIM)
ylim(YLIM)
axis square

subplot(2, 3, 6)
scatter(XYZm1(:, 1), XYZm1(:, 3), 'ko')
hold on
scatter(XYZm2(:, 1), XYZm2(:, 3), 'kx')
scatter(XYZm3(:, 1), XYZm3(:, 3), 'k+')

h = imagesc(xg, zg, reshape(Bpr./Bps, resolg, resolg)');
set(h, 'AlphaData', reshape(mask, resolg, resolg))

colormap(hot)
colorbar
xlabel('X (m)')
ylabel('Y (m)')
title('$B_p^r/ B_p^s$', 'interpreter', 'latex')

caxis([1, max(Bpr./Bps.*mask)])
xlim(XLIM)
ylim(YLIM)
axis square

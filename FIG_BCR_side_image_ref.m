clear all

Nm = 4;
Nm2 =  Nm^2;
k = 10;
Nsnaps = [100 100];

sigma2 = 0.5;

% Arrays coordinates

xx = linspace(-2, 1, Nm);

[X, Y] = meshgrid(xx, xx);

delta = 0.5;

Z1 = -2;
Z2 = 2;

XYZm1r = [X(:) Y(:)+0.5 Z1 * ones(size(X(:)))];

XYZm2r = [ Z1 * ones(size(X(:))) Y(:)+0.5 X(:)];

XYZm1 = [X(:)+delta Y(:)+0.5 Z1 * ones(size(X(:)))];

XYZm2 = [ Z1 * ones(size(X(:))) Y(:)+0.5 X(:)+delta];

% initialization grid

resolg = 51;
xg = linspace(-1.8, 1.8, resolg);
yg = 0;
zg = linspace(-1.8, 1.8, resolg);

[XG, YG, ZG] = meshgrid(xg, yg, zg);
Xg = [XG(:) YG(:) ZG(:)];

XG5 = XG(1:1:end, 1:1:end, 1:5:end);
YG5 = YG(1:1:end, 1:1:end, 1:5:end);
ZG5 = ZG(1:1:end, 1:1:end, 1:5:end);

p = 1;
%%



figure('Renderer', 'Painters');


subplot(1, 2, 1)

scatter3(XYZm1r(:, 1), XYZm1r(:, 3), XYZm1r(:, 2), 'ko')
hold on
scatter3(XYZm2r(:, 1), XYZm2r(:, 3), XYZm2r(:, 2), 'kx')

scatter3(XG5(:), ZG5(:), YG5(:), 2, [0.4, 0.4, 0.4], 'filled')
campos([5.9740  -25.6071    5.0807])

xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')

XYZs = [0.1, 0.0, 0.2];

scatter3(0.1, 0.0, 0.2, 100, [0,0,0], '+')

axis equal

subplot(1, 2, 2)

scatter3(XYZm1(:, 1), XYZm1(:, 3), XYZm1(:, 2), 'ko')
hold on
scatter3(XYZm2(:, 1), XYZm2(:, 3), XYZm2(:, 2), 'kx')

scatter3(XG5(:), ZG5(:), YG5(:), 2, [0.4, 0.4, 0.4], 'filled')
campos([5.9740  -25.6071    5.0807])

xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')

XYZs = [0.1, 0.0, 0.2];

scatter3(0.1, 0.0, 0.2, 100, [0,0,0], '+')

axis equal


%% CRBs

B3r = zeros(size(Xg, 1), 1);
B3s = zeros(size(Xg, 1), 1);
B3 = zeros(size(Xg, 1), 1);

Bps = zeros(size(Xg, 1), 1);
Bpr = zeros(size(Xg, 1), 1);
Bp = zeros(size(Xg, 1), 1);

for u = 1:size(Xg, 1)
    
[B3s(u), Bps(u), B3r(u), Bpr(u)] = BCRuncN({XYZm1, XYZm2}, Nsnaps, Xg(u, :), k,  p, sigma2, @freefieldsource);
[B3(u), Bp(u)] = BCRuncN({XYZm1r,  XYZm2r}, Nsnaps, Xg(u, :), k,  p, sigma2, @freefieldsource);

end
%%
figure('Position', [100, 100, 1000, 600])

XLIM = [-2.1 2];
YLIM = [-2.1 2];

subplot(2, 3, 1)

scatter(XYZm1r(:, 1), XYZm1r(:, 3), 'ko')
hold on
scatter(XYZm2r(:, 1), XYZm2r(:, 3), 'kx')

imagesc(xg, zg, reshape(B3, resolg, resolg)')
colormap(hot)
colorbar
xlabel('X (m)')
ylabel('Y (m)')
title('$B_{\mathbf x}^{\mathrm{ref}}$ (m$^2$)', 'interpreter', 'latex')
axis square

xlim(XLIM)
ylim(YLIM)

subplot(2, 3, 2)

scatter(XYZm1(:, 1), XYZm1(:, 3), 'ko')
hold on
scatter(XYZm2(:, 1), XYZm2(:, 3), 'kx')

imagesc(xg, zg, reshape(B3s, resolg, resolg)')
colormap(hot)
colorbar
xlabel('X (m)')
ylabel('Y (m)')
title('$B_{\mathbf x}^{\mathrm{noref}}$', 'interpreter', 'latex')

xlim(XLIM)
ylim(YLIM)
axis square

subplot(2, 3, 3)

imagesc(xg, zg, reshape(B3s./B3, resolg, resolg)');
hold on
[c, h] = contour(xg, zg, reshape(B3s./B3, resolg, resolg)', [0.6 0.8 1.0 1.2 1.4], 'k', 'ShowText', 'on', 'Color', 'w');
clabel(c, h, 'Color', 'w')

xlabel('X (m)')
ylabel('Y (m)')
title('$B_{\mathbf x}^{\mathrm{noref}}/ B_{\mathbf x}^{\mathrm{ref}}$', 'interpreter', 'latex')

axis square
axis xy


subplot(2, 3, 4)

scatter(XYZm1r(:, 1), XYZm1r(:, 3), 'ko')
hold on
scatter(XYZm2r(:, 1), XYZm2r(:, 3), 'kx')

imagesc(xg, zg, reshape(Bp, resolg, resolg)')
colormap(hot)
colorbar
xlabel('X (m)')
ylabel('Y (m)')
title('$B_p^{\mathrm{ref}}$ (Pa$^2$)', 'interpreter', 'latex')
xlim(XLIM)
ylim(YLIM)
axis square


subplot(2, 3, 5)
scatter(XYZm1(:, 1), XYZm1(:, 3), 'ko')
hold on
scatter(XYZm2(:, 1), XYZm2(:, 3), 'kx')

imagesc(xg, zg, reshape(Bps, resolg, resolg)')
colormap(hot)
colorbar
xlabel('X (m)')
ylabel('Y (m)')
title('$B_p^{\mathrm{noref}}$', 'interpreter', 'latex')
xlim(XLIM)
ylim(YLIM)
axis square

subplot(2, 3, 6)


imagesc(xg, zg, reshape(Bps./Bp, resolg, resolg)');
hold on
[cn, h] = contour(xg, zg, reshape(Bps./Bp, resolg, resolg)', [0.90:0.02:1.02], 'k', 'ShowText', 'on', 'Color', 'k');
clabel(c, h, 'Color', 'k')

xlabel('X (m)')
ylabel('Y (m)')
title('$B_p^{\mathrm{noref}}/ B_p^{\mathrm{ref}}$', 'interpreter', 'latex')

axis square
axis xy

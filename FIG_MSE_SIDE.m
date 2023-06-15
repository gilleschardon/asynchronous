clear all

figure

Nm = 4;
Nm2 = Nm^2;

K = 0.5:0.5:12;
K = 1:1:12;
NbPlots = length(K);

Nsnap = [50 200];

sigma2 = 0.5;

xx = linspace(-0.5, 0.5, Nm);
step = xx(2) - xx(1);

[X, Y] = meshgrid(xx, xx);

Z1 = -2;
Z2 = 2;

XYZm = cell(2, 1);

XYZm{1} = [X(:) Y(:) Z1 * ones(size(X(:)))];
XYZm{2} = [ Z1 * ones(size(X(:))) Y(:) X(:)];

XYZs = [-1.2, 0.0, -1.3];

LBg = [-2,-1,-2];
UBg = [0,1,0];

LB = [-2,-2,-2];
UB = [2,2,2];


% initialization grid
xg = linspace(LBg(1), UBg(1), 51);
yg = linspace(LBg(2), UBg(2), 51);
zg = linspace(LBg(3), UBg(3), 51);

[XG, YG, ZG] = meshgrid(xg, yg, zg);
Xg = [XG(:) YG(:) ZG(:)];

p = 1;

Ntest = 2000;
tic
[MSEX, MSEp] = simu_async(XYZm, XYZs, p, sigma2, Nsnap, K, LB, UB, Xg, Ntest);

for u = 1:length(K)

    k = K(u);
    [B3s(u), Bps(u), B3r(u), Bpr(u)] = BCRuncN(XYZm, Nsnap, XYZs, k,  p, sigma2, @freefieldsource);
    [B3(u), Bp(u)] = BCRuncN({cell2mat(XYZm)}, sum(Nsnap), XYZs, k,  p, sigma2, @freefieldsource);

end

save mseside1



XYZs = [-1.1, 0.0, 0.72];

LBg = [-2,-1,0];
UBg = [0,1,2];


% initialization grid
xg = linspace(LBg(1), UBg(1), 51);
yg = linspace(LBg(2), UBg(2), 51);
zg = linspace(LBg(3), UBg(3), 51);

[XG, YG, ZG] = meshgrid(xg, yg, zg);
Xg = [XG(:) YG(:) ZG(:)];


[MSEX, MSEp] = simu_async(XYZm, XYZs, p, sigma2, Nsnap, K, LB, UB, Xg, Ntest);

for u = 1:length(K)

    k = K(u);
    [B3s(u), Bps(u), B3r(u), Bpr(u)] = BCRuncN(XYZm, Nsnap, XYZs, k,  p, sigma2, @freefieldsource);
    [B3(u), Bp(u)] = BCRuncN({cell2mat(XYZm)}, sum(Nsnap), XYZs, k,  p, sigma2, @freefieldsource);

end
save mseside2
toc

%%
figure('Position', [100, 100, 800, 600])

load mseside1
plot_MSE(K, MSEX, MSEp, B3s, B3r, B3, Bps, Bpr, Bp, 2, 1, false)

load mseside2
plot_MSE(K, MSEX, MSEp, B3s, B3r, B3, Bps, Bpr, Bp, 2, 2, true)

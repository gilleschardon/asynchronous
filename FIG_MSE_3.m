clear all

figure

Nm = 4;
Nm2 =  Nm^2;

K = 1:12;
NbPlots = length(K);

Nsnap = [100, 100, 100];

sigma2 = 0.2;

xx = linspace(-0.5, 0.5, Nm);
step = xx(2) - xx(1);

[X, Y] = meshgrid(xx, xx);

Z1 = -2;
Z2 = 2;

XYZm = cell(3, 1);

XYZm{1} = [X(:) Y(:) Z1 * ones(size(X(:)))];

theta = 2*pi/3;
Rot = [cos(theta) , 0, sin(theta) ; 0, 1, 0 ; -sin(theta), 0, cos(theta)];
XYZm{2} = XYZm{1} * Rot;
XYZm{3} = XYZm{1} * Rot * Rot;

XYZs = [-1.35, 0.0, -0.2];

LBg = [-2,-1,-1];
UBg = [0,1,1];

LB = [-2,-2,-2];
UB = [2,2,2];


% initialization grid
xg = linspace(LBg(1), UBg(1), 51);
yg = linspace(LBg(2), UBg(2), 51);
zg = linspace(LBg(3), UBg(3), 51);

[XG, YG, ZG] = meshgrid(xg, yg, zg);
Xg = [XG(:) YG(:) ZG(:)];

p = 1;

Ntest = 100;

tic
[MSEX, MSEp] = simu_async(XYZm, XYZs, p, sigma2, Nsnap, K, LB, UB, Xg, Ntest);

for u = 1:length(K)

    k = K(u);
    [B3s(u), Bps(u), B3r(u), Bpr(u)] = BCRuncN(XYZm, Nsnap, XYZs, k,  p, sigma2, @freefieldsource);
    [B3(u), Bp(u)] = BCRuncN({cell2mat(XYZm)}, sum(Nsnap), XYZs, k,  p, sigma2, @freefieldsource);

end

save mse31


XYZs = [-0.43, 0.0, 0.14];

[MSEX, MSEp] = simu_async(XYZm, XYZs, p, sigma2, Nsnap, K, LB, UB, Xg, Ntest);

for u = 1:length(K)

    k = K(u);
    [B3s(u), Bps(u), B3r(u), Bpr(u)] = BCRuncN(XYZm, Nsnap, XYZs, k,  p, sigma2, @freefieldsource);
    [B3(u), Bp(u)] = BCRuncN({cell2mat(XYZm)}, sum(Nsnap), XYZs, k,  p, sigma2, @freefieldsource);

end
save mse32

toc


%%

figure('Position', [100, 100, 800, 600])

load mse31
plot_MSE(K, MSEX, MSEp, B3s, B3r, B3, Bps, Bpr, Bp, 2, 1, false)

load mse32
plot_MSE(K, MSEX, MSEp, B3s, B3r, B3, Bps, Bpr, Bp, 2, 2, true)


function plot_MSE(K, MSEX, MSEp, B3s, B3r, B3, Bps, Bpr, Bp, Nfig, nfig, leg)


MS = 10;
LW = 2;

% +xos^v*
% +++++
subplot(2, Nfig, nfig)

semilogy(K, squeeze(sum(MSEX(:, 1, :), 1)), 's', 'markersize', MS, 'linewidth', LW)
hold on
plot(K, squeeze(sum(MSEX(:, 2, :), 1)), '^', 'markersize', MS, 'linewidth', LW)
plot(K, squeeze(sum(MSEX(:, 7, :), 1)), '+', 'markersize', MS, 'linewidth', LW)
plot(K, squeeze(sum(MSEX(:, 3, :), 1)), 'x', 'markersize', MS, 'linewidth', LW)
plot(K, squeeze(sum(MSEX(:, 4, :), 1)), 'o', 'markersize', MS, 'linewidth', LW)
plot(K, squeeze(sum(MSEX(:, 6, :), 1)), 'v', 'markersize', MS, 'linewidth', LW)
plot(K, squeeze(sum(MSEX(:, 5, :), 1)), '*', 'markersize', MS, 'linewidth', LW)

plot(K, B3, '-.', 'markersize', MS, 'linewidth', LW)
plot(K, B3s, '-', 'markersize', MS, 'linewidth', LW)
plot(K, B3r, '--', 'markersize', MS, 'linewidth', LW)

if leg
    legend('MLE sync', 'MLEs', 'MLEr', 'arithm.', 'geom.', 'min', 'compl.', 'CRB sync', 'CRB strict', 'CRN relax')
end
ylabel('Position MSE (m^2)')
xlabel('k (m^{-1})')

subplot(2, Nfig, nfig+Nfig)

plot(K, MSEp(1, :), 's', 'markersize', MS, 'linewidth', LW)
hold on
plot(K, MSEp(2, :), '^', 'markersize', MS, 'linewidth', LW)
plot(K, MSEp(7, :), '+', 'markersize', MS, 'linewidth', LW)
plot(K, MSEp(3, :), 'x', 'markersize', MS, 'linewidth', LW)
plot(K, MSEp(4, :), 'o', 'markersize', MS, 'linewidth', LW)
plot(K, MSEp(6, :), 'v', 'markersize', MS, 'linewidth', LW)
plot(K, MSEp(5, :), '*', 'markersize', MS, 'linewidth', LW)

plot(K, Bp, '-.', 'markersize', MS, 'linewidth', LW)
plot(K, Bps, '-', 'markersize', MS, 'linewidth', LW)
plot(K, Bpr, '--', 'markersize', MS, 'linewidth', LW)

%legend('MLE sync', 'MLEs', 'MLEr', 'arithm.', 'geom.', 'min', 'compl.', 'B_p^c', 'B_p^s', 'B_p^r')
xlabel('k (m^{-1})')

ylabel('Power MSE (Pa^2)')



end
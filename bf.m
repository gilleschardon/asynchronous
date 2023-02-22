function bf = bf(Sigma1, x, sigma2, sv1)

% beamforming criterion (III)

% Sigma1 covariance matrix
% x position
% sigma2 noise level
% sv1 source model

g1 = sv1(x);

b1 = real(g1'*Sigma1*g1) / (sigma2 * norm(g1)^2);

bf = b1;

end

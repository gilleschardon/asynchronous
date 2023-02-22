function bf = bfprodN(Sigma, x, sigma2, sv)

% beamforming criterion, fusion by geometric averaging

% Sigma cell of covariance matrices
% x position
% sigma2 noise level
% sv cell of source models

bf = 1;

for u = 1:length(Sigma)
g = sv{u}(x);


bf = bf * real(g'*Sigma{u}*g) / (sigma2 * norm(g)^2);
end

end

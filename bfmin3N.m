function bf = bfmin3N(Sigma, S, x, sigma2, sv)

% beamforming criterion, fusion by minimum

% Sigma cell of covariance matrices
% x position
% sigma2 noise level
% sv cell of source models

bf = inf;

for u = 1:length(Sigma)

    g = sv{u}(x);

    bf = min(bf, real(g'*Sigma{u}*g) / (sigma2 * norm(g)^2));
end

end

function mlecrit = mlecritrelaxN(Sigma, S, p, x, sigma2, sv)

% MLE criterion, in function of position and powers
% for constant powers, MLE criterion of the strict model

% Sigma cell of covariance matrices
% S numbers of snapshots
% p powers
% x position
% sigma2 noise level
% sv cell of source models


mlecrit = 0;
for u = 1:length(Sigma)

g = sv{u}(x);

mlecrit = mlecrit + S(u) * (- p(u) * (real(g'*Sigma{u}*g)) / (sigma2 * (sigma2 + p(u) * norm(g)^2)) + log(sigma2 + p(u) * norm(g)^2));
end

end

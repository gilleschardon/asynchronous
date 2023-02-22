function mle = bfrelaxN(Sigma, S, x, sigma2, sv)

% MLE relaxed model, criterion in function of position

% Sigma cell of covariance matrices
% S numbers of snapshots
% x position
% sigma2 noise level
% sv cell of source models


mle = 0;

for u = 1:length(Sigma)
g = sv{u}(x);

bf = real(g'*Sigma{u}*g) /(sigma2 * norm(g)^2);

mle = mle + S(u) * (bf - log(bf));
end

end

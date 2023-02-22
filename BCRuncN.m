function [B3strict, Bpstrict, B3relax, Bprelax] = BCRuncN(Xarrays, Nsnaps, Xs, k, p, sigma2, sourcemodel)

% BCR for the strict and relaxed model

% Xarrays cell of coordinates of the microphones
% Nsnaps number of snapshots
% Xs position of the source
% k wavenumber
% p power
% sigma2 noise level
% sourcemodel sourcemodel

% B3strict BCR for position, strict model
% Bpstrict BCR for power, strict model
% B3relax BCR for position, relaxed model
% Bprelax BCR for power, relaxed model

Narrays = length(Nsnaps);

Fs = zeros(4, 4, Narrays);

for w = 1:Narrays

[a1, ax1, ay1, az1] = sourcemodel(Xarrays{w}, Xs, k);

Sigma1diff = zeros(length(a1), length(a1), 4);

% covariance matrix
Sigma1 = p * (a1 * a1') + sigma2 * eye(length(a1));

% inverse and derivatives of the covariance matrix

Sigma1inv = inv(Sigma1);
Sigma1diff(:, :, 4) = a1 * a1';
Sigma1diff(:, :, 1) = p * (a1 * ax1' + ax1 * a1');
Sigma1diff(:, :, 2) = p * (a1 * ay1' + ay1 * a1');
Sigma1diff(:, :, 3) = p * (a1 * az1' + az1 * a1');

% FIM coefficients

for u = 1:4
    for v = 1:4
        Fs(u, v, w) = Nsnaps(w) * real(trace(Sigma1inv * Sigma1diff(:, :, u) * Sigma1inv * Sigma1diff(:, :, v)));
    end
end


end

% FIM relaxed model

Frelax = zeros(3+Narrays, 3+Narrays);

for w = 1:Narrays
    Frelax(3+w, 3+w) = Fs(4, 4, w);
    Frelax(3+w, 1:3) = Fs(4, 1:3, w);
    Frelax(1:3, 3+w) = Fs(1:3, 4, w);

end

Frelax(1:3, 1:3) = sum(Fs(1:3, 1:3, :), 3);

% FIM strict model
Fstrict = sum(Fs, 3);

% BCRs

BCRrelax = inv(Frelax);
BCRstrict = inv(Fstrict);

B3strict = trace(BCRstrict(1:3, 1:3));
Bpstrict = BCRstrict(4,4);
B3relax = trace(BCRrelax(1:3, 1:3));
Bprelax = sum(sum(BCRrelax(4:end, 4:end))) / Narrays^2;



end

function Ccompl = complete_herm(Covs)

% F. Ning, J. Song, J. Hu, J. Wei, Sound source localization of non-
% synchronous measurements beamforming with block Hermitian matrix
% completion, Mechanical Systems and Signal Processing (2021)

% cannot recover the phase relationship -> do not use

UUs = cell(length(Covs), 1);

for u = 1:length(Covs)
    [U, L] = eig(Covs{u});

    UUs{u} = U * diag(sqrt(diag(L)));
end

    UU = blkdiag(UUs{:});

    MM = repmat(eye(size(Covs{1})), length(Covs), length(Covs));

    Ccompl = UU * MM * UU';

end

function [MSEX, MSEp] = simu_async(XYZm, XYZs, p, sigma2, Nsnap, K, LB, UB, Xg, Ntest)

sourcemodel = @freefieldsource;

Nasync = length(XYZm);

NbPlots = length(K);

XYZmtot = cell2mat(XYZm);

% 1 total
% 2 strict
% 3 sum
% 4 product
% 5 herm
% 6 min
% 7 relax

Xest = zeros(Ntest, 3, 7);
Pest = zeros(Ntest, 7);

MSEX = zeros(3, 7, NbPlots);
MSEp = zeros(7, NbPlots);

zzz = 0;
%% relax

for u = 1:NbPlots

    k = K(u);
    
    aasync = cell(size(XYZm));
    
    Dg = cell(size(XYZm));
    Dgnorm = cell(size(XYZm));

    dics = cell(size(XYZm));
    
    
    
    for v = 1:length(aasync)
        aasync{v} = sourcemodel(XYZm{v}, XYZs, k);
        Dg{v} =  sourcemodel(XYZm{v}, Xg, k);
        Dgnorm{v} = Dg{v} ./ sqrt(sum(abs(Dg{v}).^2, 1));
        
        dics{v} = @(x) sourcemodel(XYZm{v}, x, k);


    end
    
    
    atot = sourcemodel(XYZmtot, XYZs, k);
    Dgtot =  sourcemodel(XYZmtot, Xg, k);
    Dgnormtot = Dgtot ./ sqrt(sum(abs(Dgtot).^2, 1));


for n = 1:Ntest 
    waitbar(zzz/(NbPlots*Ntest));
    zzz = zzz+1;

    mes = cell(size(XYZm));
    Cov = cell(size(XYZm));
    
    bfgrid = zeros(size(Xg, 1), length(aasync));

    for v = 1:length(aasync)
        mes{v} = aasync{v} * p/sqrt(2) * (randn(1, Nsnap(v)) + 1i*randn(1, Nsnap(v))) + (randn(size(XYZm{v}, 1), Nsnap(v)) * sqrt(sigma2/2) + 1i*randn(size(XYZm{v}, 1), Nsnap(v)) * sqrt(sigma2/2));
        
        Cov{v} = mes{v} * mes{v}' / Nsnap(v);
        
        bfgrid(:, v) = real(sum((Dgnorm{v}' * Cov{v}) .* Dgnorm{v}.', 2)) / sigma2;
    end
    mestot = atot * p/sqrt(2) * (randn(1, sum(Nsnap)) + 1i*randn(1, sum(Nsnap))) + (randn(size(XYZmtot, 1), sum(Nsnap)) * sqrt(sigma2/2) + 1i*randn(size(XYZmtot, 1), sum(Nsnap)) * sqrt(sigma2/2));
    Covtot = mestot * mestot' / sum(Nsnap);


    % sum
    
    sumbf = sum(bfgrid, 2);
    [~, idxsum] = max(sumbf);
    

    objsum = @(x) - bfsumN(Cov, x, sigma2, dics);
    Xsum = fmincon(@(x)objsum(x), Xg(idxsum, :), [], [], [], [], LB, UB);
    
    Xest(n, :, 3) = Xsum;
    
    Pr = zeros(Nasync, 1);
    
    for v = 1:Nasync   
       gr = dics{v}(Xsum);

        Pr(v) = (real(gr'*Cov{v}*gr)/norm(gr)^2 - sigma2) / norm(gr)^2;
    end
    
    
    Pest(n, 3) = mean(Pr);
    
    % MLE relax
    
    mlerbf = sum(bfgrid - log(bfgrid), 2);
    [~, idxmler] = max(mlerbf);
    

    objmler = @(x) - bfrelaxN(Cov, Nsnap, x, sigma2, dics);
    Xmler = fmincon(@(x)objmler(x), Xg(idxmler, :), [], [], [], [], LB, UB);
    
    Xest(n, :, 7) = Xmler;
    
    Pmler = zeros(Nasync, 1);
    
    for v = 1:Nasync   
       gr = dics{v}(Xmler);

        Pmler(v) = (real(gr'*Cov{v}*gr)/norm(gr)^2 - sigma2) / norm(gr)^2;
    end
    
    
    Pest(n, 7) = mean(Pmler);

    
    
    
    % MLE strict
    
    % the strict objective is the relaxed objective with same power for
    % each array
    objstrict = @(x) mlecritrelaxN(Cov, Nsnap, repmat(x(4), Nasync, 1), x(1:3), sigma2, dics);
    XPstrict = fmincon(@(x)objstrict(x), [Xsum mean(Pr)], [], [], [], [], [LB 0], [UB inf]);
    
    Xest(n, :, 2) = XPstrict(1:3);
    
    Pest(n, 2) = XPstrict(4);

    % prod
    
    prodbf = prod(bfgrid, 2);
    [~, idxprod] = max(prodbf);
    
    objprod = @(x) - bfprodN(Cov, x, sigma2, dics);
    Xprod = fmincon(@(x)objprod(x), Xg(idxprod, :), [], [], [], [], LB, UB);
    
    Xest(n, :, 4) = Xprod;
    
    for v = 1:Nasync   
       gr = dics{v}(Xsum);

        Pps(v) = (real(gr'*Cov{v}*gr)/norm(gr)^2 - sigma2) / norm(gr)^2;
    end
  
  
     Pest(n, 4) = (prod(Pps))^(1/Nasync);

     % hermi
     
    Covherm = complete_herm(Cov);
 
     bfherm = real(sum((Dgnormtot' * Covherm) .* Dgnormtot.', 2));
     [~, idxherm] = max(bfherm);
     objbf = @(x) - bf(Covherm, x, sigma2, @(x) sourcemodel(XYZmtot, x, k));
 
     Xherm = fmincon(@(x)objbf(x), Xg(idxherm, :), [], [], [], [], LB, UB);
     Xest(n, :, 5) = Xherm;
     
     gh = sourcemodel(XYZmtot, Xherm, k);
 
     Pest(n, 5) = (real(gh'*Covherm*gh)/norm(gh)^2 - sigma2) / norm(gh)^2;
    
      

    
     % min3
     
    minbf3g = min(bfgrid, [], 2);
    [~, idxminbf3] = max(minbf3g);
    
    objmin = @(x) - bfmin3N(Cov, Nsnap, x, sigma2, dics);
    Xmin  = fmincon(@(x)objmin(x), Xg(idxminbf3, :), [], [], [], [], LB, UB);
    
    Xest(n, :, 6) = Xmin;

        
    Pms = zeros(Nasync, 1);
    for v = 1:Nasync   
       gm = dics{v}(Xmin);

        Pms(v) = (real(gm'*Cov{v}*gm)/norm(gm)^2 - sigma2) / norm(gm)^2;
    end
    Pest(n, 6) = min(Pms);
    % total

    bfgridtot = real(sum((Dgnormtot' * Covtot) .* Dgnormtot.', 2));
    [~, idxtot] = max(bfgridtot);
    objbf = @(x) - bf(Covtot, x, sigma2, @(x) sourcemodel(XYZmtot, x, k));

    Xtot = fmincon(@(x)objbf(x), Xg(idxtot, :), [], [], [], [], LB, UB);
    
    Xest(n, :, 1) = Xtot;
    
    gh = sourcemodel(XYZmtot, Xtot, k);
    Pest(n, 1) = (real(gh'*Covtot*gh)/norm(gh)^2 - sigma2) / norm(gh)^2;

end

MSEX(:, :, u) = squeeze(mean((Xest - XYZs).^2, 1));
MSEp(:, u) = squeeze(mean((Pest - p).^2, 1));

end


end
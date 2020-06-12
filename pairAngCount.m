function [N_corr,N_all] = pairAngCount(X,th_ed,verbose)
% PAIRANGCOUNT counts difference angles between pairs.
% DKS 2020

if nargin<3
    verbose = 0;
end


%%% setup data
nShot = size(X,1);
M = size(X,2);

X1 = X(:,1);
if M > 2
    error('X must be Nx1 or Nx2 cell araray.');
elseif M == 1
    X2 = X1;
else
    X2 = X(:,2);
end


%%% pair counting
nbin_th = length(th_ed)-1;
N_corr = zeros(nbin_th,1);
N_all = N_corr;

if verbose>0
    progressbar('pairAngCount');
end
for ii = 1:nShot
    for jj = 1:nShot
        %%% pairwise diff angle
        if (ii == jj) && (M == 1)
            % self-counting
            dtheta = pair_diff_angle(X1{ii});
        else
            dtheta = pair_diff_angle(X1{ii},X2{jj});
        end
        
        %%% histogram
        tN = nhist(dtheta(:),{th_ed});
        % update
        if ii == jj
            N_corr = N_corr + tN;
        end
        N_all = N_all + tN;
    end
    if verbose>0
        progressbar(ii/nShot);
    end
end

end
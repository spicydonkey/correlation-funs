function [N_dth_corr,N_dth_all] = countPairDiffAngles(X,th_ed,verbose)
% countPairDiffAngles counts difference angles between pairs.
%
%   [N_corr,N_all] = countPairDiffAngles(X,th_ed,verbose)
%
% X             Experimental data.  Cell-array (n_shot x m) of shots (n_counts x 3 array).
% th_ed         Bin edges (1D vector) to histogram difference angles.  
% 
% N_dth_corr    Total number of pairwise difference angles from same shot (correlated).
% N_dth_all     Total number of pairwise difference angles from all data.
%
% DKS 2020

if nargin<3
    verbose = 1;
end

%%% setup data
n_shot = size(X,1);
M = size(X,2);

X1 = X(:,1);
if M > 2
    error('X must be Nx1 or Nx2 cell araray.');
elseif M == 1
    X2 = X1;
else
    X2 = X(:,2);
end

%% pair counting
nbin_th = length(th_ed)-1;
N_dth_corr = zeros(nbin_th,1);
N_dth_all = N_dth_corr;

if verbose>0
    progressbar('pairAngCount');
end

if M == 1
    %% single species -- symmetry
    counter = 0;
    n_pairshot_cases = n_shot*(n_shot+1)/2;       % # unordered pairs {ii,jj} | ii,jj = 1...nShot
    for ii = 1:n_shot
        for jj = ii:n_shot       % symmetry between pair of shots (ii,jj) = (jj,ii)
            % get angles between every pair (except self-pairing)
            if ii == jj
                dtheta = pair_diff_angle(X{ii});
            else
                dtheta = pair_diff_angle(X{ii},X{jj});
            end
            
            % histogram difference angles
            this_N = nhist(dtheta(:),{th_ed});
            if ii == jj
                % pairs observed same shot
                N_dth_corr = N_dth_corr + this_N;
                N_dth_all = N_dth_all + this_N;
            else
                % pairs observed in different shots
                N_dth_all = N_dth_all + this_N;       
            end
            
            counter = counter + 1;
            if verbose>0
                progressbar(counter/n_pairshot_cases);
            end
        end
    end
else
    %% two species: no symmetry (ii,jj) ~= (jj,ii)
    for ii = 1:n_shot
        for jj = 1:n_shot
            dtheta = pair_diff_angle(X1{ii},X2{jj});
            
            this_N = nhist(dtheta(:),{th_ed});
            if ii == jj
                N_dth_corr = N_dth_corr + this_N;
            end
            N_dth_all = N_dth_all + this_N;
        end
        if verbose>0
            progressbar(ii/n_shot);
        end
    end
end
end
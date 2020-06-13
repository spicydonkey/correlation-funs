function [g2,N_dth_corr,N_dth_all] = angular_g2(K,th_ed,verbose)
%ANGULAR_G2 evaluates g2 correlation function for relative difference angle
%between atoms.
%
% X             Experimental data.  1D cell-array (n_shot x 1) of shots (n_counts x 3 array).
% th_ed         Bin edges (1D vector) to histogram difference angles.  
% 
% g2            First particle-integrated g2 function of relative angle.
% N_dth_corr    Total number of pairwise difference angles from same shot (correlated).
% N_dth_all     Total number of pairwise difference angles from all data.
% 
if nargin<3
    verbose = 1;
end
K = unpackrows(K);
[N_dth_corr,N_dth_all] = countPairDiffAngles(K,th_ed,verbose);
% normalise
nn_dth_corr = N_dth_corr./sum(N_dth_corr);
nn_dth_all = N_dth_all./sum(N_dth_all);
g2 = nn_dth_corr./nn_dth_all;
end
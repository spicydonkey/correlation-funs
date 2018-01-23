function [n,N] = ncorr_ang(ks,dth_ed)
%Boilerplate 2-particle normalised coincidence counter
%
%   ks: Nx3 array of sph-polar vecs
%   dth_ed: (nbin+1)x1 array; diff-theta bin edges
%
%   n: rel-angle coincidence rates, normalised to unit-sum
%   N: hist-counts of coincidences

nCounts=size(ks,1);
nbins=numel(dth_ed)-1;
N=zeros(nbins);

% histogram the rel-angles between "all" unique atom pairs
for ii=1:nCounts
    this_dth=diffAngleSph(ks(ii,1),ks(ii,2),...
            ks((ii+1):end,1),ks((ii+1):end,2));
    N=N+nhist(this_dth,{dth_ed});
end
% normalise histogram
nTotPairs=nCk(nCounts,2);
n=N/nTotPairs;      % normalised to unit-sum

end
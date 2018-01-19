function [g2,G2_corr,G2_unc]=g2_ang(k,dth_ed)
%2-particle correlations in relative angle in spherical-polar coord
%
%   k: nShotx1 cell-array of cart ZXY coords
%   dth_ed: (nbin+1)x1 array; diff-theta bin edges
%
%   g2: 1D angular correlation
%   G2_corr: unnormalised coincidences in same experiment
%   G2_unc: like above but from different shots
%
%
%   NOTE:
%       correlation is ONLY evaluted in difference angle.
%       vector norms are ignored.
%

% cart --> sph-polar
ks=cellfun(@(x) zxy2sphpol(x),k,'UniformOutput',false);

nShot=size(ks,1);
nCounts=shotSize(ks);
nbins=numel(dth_ed)-1;

%%% correlated 2-particle histogram
G2_corr=zeros(nbins,1);
parfor ii=1:nShot
    this_ks=ks{ii};     % this shot
    this_n=nCounts(ii); % number of atoms in this shot
    
    for jj=1:this_n
        % rel-angle between all unique atom-pairs
        this_dth=diffAngleSph(this_ks(jj,1),this_ks(jj,2),...
            this_ks((jj+1):end,1),this_ks((jj+1):end,2));
        G2_corr=G2_corr+nhist(this_dth,{dth_ed});     % update G2 hist
    end
end

%%% uncorrelated 2-particle histogram
% TODO - this scales horribly with number of shots
G2_unc=zeros(nbins,1);
parfor ii=1:nShot
    % self-included in "uncorrelated"
    this_ks_collated=vertcat(ks{ii:end});       % TODO - this makes entire data a broadcast variable for parfor
    this_ks=ks{ii};     % ref parts
    this_n=nCounts(ii);
    
    for jj=1:this_n
        this_dth=diffAngleSph(this_ks(jj,1),this_ks(jj,2),...
            this_ks_collated((jj+1):end,1),this_ks_collated((jj+1):end,2));
        G2_unc=G2_unc+nhist(this_dth,{dth_ed});
    end
end

%%% normalise correlation function
% pair combinatorics
%   count all unique pairs from data
nPairsCorr=sum(arrayfun(@(n)nCk(n,2),nCounts));
nPairsUnc=nCk(sum(nCounts),2);

% normalise G2
G2_corr=G2_corr/nPairsCorr;
G2_unc=G2_unc/nPairsUnc;

% NOTE: pre-filtering can significantly improve SNR
g2=G2_corr./G2_unc;
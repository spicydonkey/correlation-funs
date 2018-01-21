function [g2,G2_corr,G2_unc]=g2x_ang(k,dth_ed)
%Two-species g(2) analysis in angular coords
%
%   k: nShotx2 cell-array of cart ZXY coords
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
    this_ks=ks{ii,1};       % first atoms
    this_qs=ks{ii,2};       % second atoms
    this_n=nCounts(ii,1);   % number of first atoms (refs)
    
    for jj=1:this_n
        % X-species, rel-angle between all unique atom-pairs
        this_dth=diffAngleSph(this_ks(jj,1),this_ks(jj,2),...
            this_qs(:,1),this_qs(:,2));
        G2_corr=G2_corr+nhist(this_dth,{dth_ed});     % update G2 hist
    end
end

%%% uncorrelated 2-particle histogram
% TODO - this scales horribly with number of shots
G2_unc=zeros(nbins,1);
qs_collated=vertcat(ks{:,2});   % all second species collated
parfor ii=1:nShot
    this_ks=ks{ii,1};       % first particle
    this_n=nCounts(ii,1);
    
    for jj=1:this_n
        this_dth=diffAngleSph(this_ks(jj,1),this_ks(jj,2),...
            qs_collated(:,1),qs_collated(:,2));
        G2_unc=G2_unc+nhist(this_dth,{dth_ed});
    end
end

%%% normalise correlation function
% x-species pair combinatorics
nPairsCorr=sum(prod(nCounts,2));    % no. unique pairs x-species same exp
nPairsUnc=prod(sum(nCounts,1));     % no. unique pairs x-species all exp

% normalise G2
G2_corr=G2_corr/nPairsCorr;
G2_unc=G2_unc/nPairsUnc;

% NOTE: pre-filtering can significantly improve SNR
g2=G2_corr./G2_unc;

end
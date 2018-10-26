function [g2,G2_shot,G2_norm]=g2x_ang(k,dth_ed)
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
G2_shot=zeros(nbins,1);
parfor ii=1:nShot
    tks=ks{ii,1};       % first atoms
    tqs=ks{ii,2};       % second atoms
    tn=nCounts(ii,1);   % number of first atoms (refs)
    
    for jj=1:tn
        % X-species, rel-angle between all unique atom-pairs
        tdth=diffAngleSph(tks(jj,1),tks(jj,2),...
            tqs(:,1),tqs(:,2));
        G2_shot=G2_shot+nhist(tdth,{dth_ed});     % update G2 hist
    end
end

%%% uncorrelated 2-particle histogram
% TODO - this scales horribly with number of shots
G2_norm=zeros(nbins,1);
qs_collated=vertcat(ks{:,2});   % all second species collated
parfor ii=1:nShot
    tks=ks{ii,1};       % first particle
    tn=nCounts(ii,1);
    
    for jj=1:tn
        tdth=diffAngleSph(tks(jj,1),tks(jj,2),...
            qs_collated(:,1),qs_collated(:,2));
        G2_norm=G2_norm+nhist(tdth,{dth_ed});
    end
end

%%% normalise correlation function
% x-species pair combinatorics
nPairsShot=sum(prod(nCounts,2));    % no. unique pairs x-species same exp
nPairsNorm=prod(sum(nCounts,1));     % no. unique pairs x-species all exp

% normalise G2
G2_shot=G2_shot/nPairsShot;
G2_norm=G2_norm/nPairsNorm;

% NOTE: pre-filtering can significantly improve SNR
g2=G2_shot./G2_norm;

end
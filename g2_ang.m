function [g2,G2_shot,G2_norm]=g2_ang(k,dth_ed)
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
G2_shot=zeros(nbins,1);
parfor ii=1:nShot
    tks=ks{ii};     % this shot
    tn=nCounts(ii); % number of atoms in this shot
    for jj=1:tn
        % rel-angle between all unique atom-pairs
        tdth=diffAngleSph(tks(jj,1),tks(jj,2),...
            tks((jj+1):end,1),tks((jj+1):end,2));
        G2_shot=G2_shot+nhist(tdth,{dth_ed});     % update G2 hist
    end

%     % encapsulated the g2 boilerplate
%     [~,Ncorr]=ncorr_ang(this_ks,dth_ed);
%     G2_corr=G2_corr+Ncorr;
end

%%% uncorrelated 2-particle histogram
% TODO - this scales horribly with number of shots
G2_norm=zeros(nbins,1);
parfor ii=1:nShot
    % self-included in "uncorrelated"
    this_ks_collated=vertcat(ks{ii:end});       % TODO - this makes entire data a broadcast variable for parfor
    tks=ks{ii};     % ref parts
    tn=nCounts(ii);
    
    for jj=1:tn
        tdth=diffAngleSph(tks(jj,1),tks(jj,2),...
            this_ks_collated((jj+1):end,1),this_ks_collated((jj+1):end,2));
        G2_norm=G2_norm+nhist(tdth,{dth_ed});
    end
end

% % encapsulated the g2 boilerplate
% [~,Nunc]=ncorr_ang(vertcat(ks{:}),dth_ed);
% % NB the above is going to be computationally expensive: uncorrelated G2 
% %   by collating all shots and running 2-part coincidence histogram
% G2_unc=Nunc;

%%% normalise correlation function
% pair combinatorics
%   count all unique pairs from data
nPairsShot=sum(arrayfun(@(n)nCk(n,2),nCounts));
nPairsNorm=nCk(sum(nCounts),2);

% normalise G2
G2_shot=G2_shot/nPairsShot;
G2_norm=G2_norm/nPairsNorm;

% NOTE: pre-filtering can significantly improve SNR
g2=G2_shot./G2_norm;

end
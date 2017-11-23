function [g2,G2_shot,G2_norm]=g2x_bb(k,dk_ed)
% Two-species g(2) analysis in Cartesian coords
% * particle correlations *across* species
%
% [g2, G2_shot, G2_norm] = g2x_bb(k, dk_ed)
%
% 
% k: nShot x 2 cell array of particle positions (cart zxy)
% dk_ed: diff-k bin edges (1x3 Cart cell-array of (ncent+1)x1 array)
%
% NOTE
% * ensure function args format. no error checking is implemented
% * see g2_bb.m for single-species g2
%
% TODO
%

nShot=size(k,1);
nCounts=cellfun(@(x) size(x,1),k);

ndk_bins=cellfun(@(ed) numel(ed)-1,dk_ed);


%%% shot-to-shot 2-particle histogram
G2_shot=zeros(ndk_bins);
parfor ii=1:nShot
    tk=k{ii,1};         % this shot
    tq=k{ii,2};         % target to search pairs
    tn=nCounts(ii,1);   % number of atoms in this shot
    
    for jj=1:tn
        % X-species BB condition
        ttK=tk(jj,:);   % this reference atom
        tdK=tq+ttK;     % BB diff-vecs to all unique paris
        G2_shot=G2_shot+nhist(tdK,dk_ed);   % update G2 hist
    end
end

%%% normalising 2-particle histogram
% TODO - scaling problem - can optimise for speed
G2_norm=zeros(ndk_bins);
tq_collated=vertcat(k{:,2});    % all 2nd (target) species collated
parfor ii=1:nShot
    tk=k{ii,1};     % this ref shot
    
    % NOTE boilerplate
    tn=nCounts(ii,1);
    
    for jj=1:tn
        % X-species BB condition
        ttK=tk(jj,:);
        tdK=tq_collated+ttK;
        G2_norm=G2_norm+nhist(tdK,dk_ed);
    end
end

%%% normalise correlation function
% x-species pair combinatorics
nPairsShot=sum(prod(nCounts,2));    % no. unique pairs x-species same exp
nPairsNorm=prod(sum(nCounts,1));    % no. unique pairs x-species all exp


% normalise G2
G2_shot=G2_shot/nPairsShot;
G2_norm=G2_norm/nPairsNorm;

% TODO - can smooth before normalising
g2=G2_shot./G2_norm;

end

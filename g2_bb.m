function [g2,G2_shot,G2_norm]=g2_bb(k,dk_ed)
% Single component g(2) analysis in Cartesian coords
%
% [g2, G2_shot, G2_norm] = g2_bb(k, dk_ed)
%
% * this is the simplest way to call back-to-back g2 correlation
%
%
% k: nShot x 1 cell array of *halo* counts
% dk_ed: diff-k bin edges (1x3 Cart cell-array of (ncent+1)x1 array)
%
% TODO
% * [ ] boilerplate code in G2
%   

nShot=size(k,1);

ndk_bins=cellfun(@(ed)numel(ed)-1,dk_ed);


%%% shot-to-shot 2-particle histogram
G2_shot=zeros(ndk_bins);
parfor ii=1:nShot
    tk=k{ii};       % this shot
    tn=size(tk,1);  % number of atoms in this shot
    
%     for jj=1:(tn-1)
    for jj=1:tn         % TODO - check if this OK
        % BB condition
        ttK=tk(jj,:);   % this atom which we are looking BB partners for
        tdK=tk((jj+1):end,:)+ttK;   % BB diff-vectors to unique pairs
        G2_shot=G2_shot+nhist(tdK,dk_ed);   % update G2 hist
    end
end

%%% normalising 2-particle histogram
% TODO - this scales horribly with number of shots
G2_norm=zeros(ndk_bins);
parfor ii=1:nShot
    tk_collated=vertcat(k{ii:end});     % TODO - this makes entire `k` a broadcast variable for parfor
    
    % NOTE here is the boilerplate
    tn=size(tk_collated,1);
    
    for jj=1:tn
        % BB condition
        ttK=tk_collated(jj,:);   % this atom which we are looking BB partners for
        tdK=tk_collated((jj+1):end,:)+ttK;   % BB diff-vectors to unique pairs
        G2_norm=G2_norm+nhist(tdK,dk_ed);   % update G2 hist
    end
end

%%% normalise correlation function
nCounts=cellfun(@(x) size(x,1),k);      % number of atoms in each shot

% pair combinatorics
% TODO Check maths
nPairsShot=sum(nCounts.*(nCounts-1)/2);     % no. unique pairs in shot
nPairsNorm=(sum(nCounts)^2-sum(nCounts.^2))/2;  % no. unique pairs in dataset

% normalise G2
G2_shot=G2_shot/nPairsShot;
G2_norm=G2_norm/nPairsNorm;

% TODO - can smooth before normalising
g2=G2_shot./G2_norm;
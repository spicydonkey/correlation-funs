function nn = ncorrx_bb(k1,k2,dk_ed)
% Boilerplate normalised back-to-back (BB) pair x-correlation counter
%
% nn = ncorr_bb(k1,k2,dk_ed)
%
% k1,k2:    N x 3 array of k-vectors
% dk_ed:    diff-k bin edges (1x3 Cart cell-array
%

nCounts1=size(k1,1);     % number of atoms
nCounts2=size(k2,1);     % number of atoms
ndk_bins=cellfun(@(ed)numel(ed)-1,dk_ed);

nn=zeros(ndk_bins);
% X-state BB condition for all unique pairs
for ii=1:nCounts1
    tK=k1(ii,:);        % this atom which we are looking x-state BB partners for
    dK=k2+tK;           % x-state BB diff-vectors to unique pairs
    nn=nn+nhist(dK,dk_ed);      % update correlation histogram
end
% normalise histogram
nTotPairs=nCounts1*nCounts2;    % no. unique pairs x-species
nn=nn/nTotPairs;

end
function nn = ncorr_bb(k,dk_ed)
% Boilerplate normalised back-to-back (BB) pair correlation counter
%
% nn = ncorr_bb(k,dk_ed)
%
% k:        N x 3 array of k-vectors
% dk_ed:    diff-k bin edges (1x3 Cart cell-array
%

nCounts=size(k,1);     % number of atoms
ndk_bins=cellfun(@(ed)numel(ed)-1,dk_ed);

nn=zeros(ndk_bins);
% BB condition for all unique pairs
for ii=1:nCounts
    tK=k(ii,:);     % this atom which we are looking BB partners for
    dK=k((ii+1):end,:)+tK;      % BB diff-vectors to unique pairs
    nn=nn+nhist(dK,dk_ed);      % update correlation histogram
end
% normalise histogram
nTotPairs=nCk(nCounts,2);
nn=nn/nTotPairs;

end
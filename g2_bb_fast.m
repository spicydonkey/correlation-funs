function [g2,G2_shot,G2_norm,G2_shot_sdev,G2_norm_sdev]=g2_bb_fast(k,dk_ed,rsamp)
% Fast Single component g(2) analysis in Cartesian coords
%   correlations from uncorrelated random events are truncated by random sampling from all shots
%
% [g2, G2_shot, G2_norm] = g2_bb_fast(k, dk_ed, rsamp)
%
% * this is the simplest way to call back-to-back g2 correlation
%
%
% k:        nShot x 1 cell array of particle positions (cart zxy)
% dk_ed:    diff-k bin edges (1x3 Cart cell-array of (ncent+1)x1 array)
% rsamp:    fraction of total uncorrelated data to sample for speed-up in normalising g2
%
% TODO
% * [ ] return G2 grid centers

nShot=size(k,1);
nCounts=shotSize(k);
ndk_bins=cellfun(@(ed)numel(ed)-1,dk_ed);

%%% shot-to-shot 2-particle histogram
ncorr_shot=zeros([ndk_bins,nShot]);     % preallocate correlation matrix
parfor ii=1:nShot
    tk=k{ii};       % k-vecs in this shot
    ncorr_shot(:,:,:,ii)=ncorr_bb(tk,dk_ed);        % update G2 hist
end
G2_shot=mean(ncorr_shot,4,'omitnan');           % shot averaged G2 (omitnan is required to deal with empty halos)
G2_shot_sdev=std(ncorr_shot,0,4,'omitnan');

%%% normalising 2-particle histogram
k_coll=vertcat(k{1:end});       % all shots collated
nTotCounts=sum(nCounts);
nSampCounts=round(nTotCounts*rsamp);

nRepSamp=round(1/rsamp);                % no. time to repeat fast uncorrelated G2
ncorr_norm=zeros([ndk_bins,nRepSamp]);  % preallocate correlation matrix
parfor ii=1:nRepSamp
    % TODO - does the randperm truly generate random samples when parallelised?
    tiSamp=randperm(nTotCounts,nSampCounts);     % randomly sample from all counts
    tk_coll_samp=k_coll(tiSamp,:);        % get the sampled set of k-vecs
    
    ncorr_norm(:,:,:,ii)=ncorr_bb(tk_coll_samp,dk_ed);   % update G2 hist
end
G2_norm=mean(ncorr_norm,4,'omitnan');       % rep averaged G2 (omitnan is required to deal with empty halos)
G2_norm_sdev=std(ncorr_norm,0,4,'omitnan');


% TODO - can smooth before normalising
g2=G2_shot./G2_norm;

end
function [g2,G2_shot,G2_norm,G2_shot_sdev,G2_norm_sdev]=g2x_bb_fast(k,dk_ed,rsamp)
% Fast two-species g(2) analysis in Cartesian coords
% * particle correlations *across* species
% * correlations from uncorrelated random events are truncated by random sampling from all shots
%
% [g2,G2_shot,G2_norm,G2_shot_sdev,G2_norm_sdev]=g2x_bb_fast(k,dk_ed,rsamp)
%
% 
% k: nShot x 2 cell array of particle positions (cart zxy)
% dk_ed: diff-k bin edges (1x3 Cart cell-array of (ncent+1)x1 array)
%
% NOTE
% * ensure function args format. no error checking is implemented
% * see g2_bb_fast.m for single-species g2
%
% TODO
% * [ ] return G2 grid centers

nShot=size(k,1);
nCounts=shotSize(k);
ndk_bins=cellfun(@(ed) numel(ed)-1,dk_ed);

%%% shot-to-shot 2-particle histogram
ncorr_shot=zeros([ndk_bins,nShot]);     % preallocate correlation matrix
parfor ii=1:nShot
    tk=k{ii,1};         % this shot
    tq=k{ii,2};         % target to search pairs
    
    ncorr_shot(:,:,:,ii)=ncorrx_bb(tk,tq,dk_ed);        % update G2 hist
end
G2_shot=mean(ncorr_shot,4,'omitnan');           % shot averaged G2 (omitnan is required to deal with empty halos)
G2_shot_sdev=std(ncorr_shot,0,4,'omitnan');


%%% normalising 2-particle histogram
k_coll=vertcat(k{:,1});                     % counts collated 
q_coll=vertcat(k{:,2});    
nTotCounts1=sum(nCounts(:,1));              % total counts across shots
nTotCounts2=sum(nCounts(:,2));
nSampCounts1=round(nTotCounts1*rsamp);      % no. counts to sample for decimated g2
nSampCounts2=round(nTotCounts2*rsamp);

nRepSamp=round(1/rsamp);                % no. time to repeat fast uncorrelated G2
ncorr_norm=zeros([ndk_bins,nRepSamp]);  % preallocate correlation matrix
parfor ii=1:nRepSamp
    % TODO - does the randperm truly generate random samples when parallelised?
    tiSamp1=randperm(nTotCounts1,nSampCounts1);     % randomly sample from all counts
    tiSamp2=randperm(nTotCounts2,nSampCounts2);     % randomly sample from all counts
    
    tk_coll_samp=k_coll(tiSamp1,:);        % get the sampled set of k-vecs
    t1_coll_samp=q_coll(tiSamp2,:);        % get the sampled set of q-vecs
    
    ncorr_norm(:,:,:,ii)=ncorrx_bb(tk_coll_samp,t1_coll_samp,dk_ed);   % update G2 hist
end
G2_norm=mean(ncorr_norm,4,'omitnan');       % rep averaged G2 (omitnan is required to deal with empty halos)
G2_norm_sdev=std(ncorr_norm,0,4,'omitnan');

% TODO - can smooth before normalising
g2=G2_shot./G2_norm;

end

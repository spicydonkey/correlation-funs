%%% create halos
nShots=200;
nPairs=10;
k_dither=[0.03,0.03,0.03];
det_qe=1;

k_dist=cell(nShots,2);
for ii=1:nShots
    [k_dist{ii,1},k_dist{ii,2}]=genCorrSource(nPairs,k_dither,det_qe);
end

% partners combined and indistinguishable
k=cellfun(@(k1,k2) vertcat(k1,k2),k_dist(:,1),k_dist(:,2),'UniformOutput',false);

%%% TEST g2
% set up bins
dk_ed_vec=linspace(-0.3,0.3,20);
dk_cent_vec=dk_ed_vec(1:end-1)+0.5*diff(dk_ed_vec);
dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};
dk_cent={dk_cent_vec,dk_cent_vec,dk_cent_vec};

dk=ndgrid(dk_cent{:});      % grid of dk centers

% run g2
[g2,G2s,G2n]=g2_bb(k,dk_ed);

% plot g2
h=figure(1);
imagesc(dk_cent_vec,dk_cent_vec,g2(:,:,10));
colorbar();
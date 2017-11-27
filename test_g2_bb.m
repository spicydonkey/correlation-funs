%%% create halos
nShots=1e4;
nPairs=100;
k_dither=0.0141*[1,1,1];
det_qe=0.1;

k_src=cell(nShots,2);
for ii=1:nShots
    [k_src{ii,1},k_src{ii,2}]=genCorrSource(nPairs,k_dither,det_qe);
end

% partners combined and indistinguishable
k=cellfun(@(k1,k2) vertcat(k1,k2),k_src(:,1),k_src(:,2),'UniformOutput',false);

%% TEST g2
% set up bins
nbin=30;        % keep it even
dk_ed_vec=linspace(-0.2,0.2,nbin);      % keep it symmetric
dk_cent_vec=dk_ed_vec(1:end-1)+0.5*diff(dk_ed_vec);
dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};  % keep it symmetric
dk_cent={dk_cent_vec,dk_cent_vec,dk_cent_vec};

dk=ndgrid(dk_cent{:});      % grid of dk centers

% run g2
% [g2,G2s,G2n]=g2_bb(k,dk_ed);          % simple BB particle counts
% [g2,G2s,G2n]=g2_bb_fast(k,dk_ed,0.1);          % simple BB particle counts

% [g2,G2s,G2n]=g2x_bb(k_src,dk_ed);      % 2-species BB
[g2,G2s,G2n]=g2x_bb_fast(k_src,dk_ed,0.1);      % 2-species BB

%% plot g2
idx_zero=ceil(nbin/2);

h_2d=figure();
imagesc(dk_cent_vec,dk_cent_vec,g2(:,:,idx_zero));
cbar=colorbar();
cbar.Label.String='g^{(2)}_{BB}';
xlabel('\Delta k_{i}');
ylabel('\Delta k_{j}');


h1d=figure();
hold on;
plot(dk_cent_vec,squeeze(g2(:,idx_zero,idx_zero)),'o-');
plot(dk_cent_vec,squeeze(g2(idx_zero,:,idx_zero)),'^-');
plot(dk_cent_vec,squeeze(g2(idx_zero,idx_zero,:)),'*-');
xlabel('\Delta k');
ylabel('g^{(2)}_{BB}');

box on;
grid on;

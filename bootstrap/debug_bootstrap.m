%% DEBUG: test from here
n_frac_samp=0.05;
n_subset=160;                    % no. of bootstrap repeats


%% PRELIM bootstrapping
%   TODO
%   [x] with replacement
%   [ ] convergence
%

%%% CONFIG
% DEBUG (2 lines below commented out)
% n_frac_samp=1/5;
% n_subset=10;                    % no. of bootstrap repeats
%   NOTE: unclear at the moment how config affects analysis

% dataset and subset
nshot_par=shotSize(k_par);    % num. exp shots for each scanned parameter set
subset_shotsize=round(nshot_par*n_frac_samp);    % shot-size of bootstrap sampled subset

g2_bootstrap=cell(nparam,1);
g2anti_samp=cell(nparam,1);
g2corr_samp=cell(nparam,1);
E_samp=cell(nparam,1);
E0_samp=cell(nparam,1);

for ii=1:nparam
    tk_par=k_par{ii};   % the population-representative data set
    
    % random sub-sample sel with replacement
    Isamp=randi(nshot_par(ii),[n_subset,subset_shotsize(ii)]);
    
    % initialise
    g2_bootstrap{ii}=cell(1,3);
    for jj=1:3
        g2_bootstrap{ii}{jj}=NaN([dk_grid_size,n_subset]);
    end
    g2anti_samp{ii}=NaN(n_subset,1);
    g2corr_samp{ii}=NaN(n_subset,1);
    E_samp{ii}=NaN(n_subset,1);
    E0_samp{ii}=NaN(n_subset,1);
    
    for jj=1:n_subset
        k_samp=tk_par(Isamp(jj,:),:);           % the sampled data
        
        tg2=summary_disthalo_g2(k_samp,dk_ed,0,0,0,0);      % evaluate function
        
        for kk=1:3
            g2_bootstrap{ii}{kk}(:,:,:,jj)=tg2{kk};
        end
        
        % *aproximate* g2 amplitude at evaluated value at dk=0
        tg2corr=mean([tg2{1}(idx_dk0,idx_dk0,idx_dk0),tg2{2}(idx_dk0,idx_dk0,idx_dk0)]);     % get results
        tg2anti=tg2{3}(idx_dk0,idx_dk0,idx_dk0);
        
        [tE,tE0]=g2toE(tg2corr,tg2anti);
        
        % store
        g2anti_samp{ii}(jj)=tg2anti;
        g2corr_samp{ii}(jj)=tg2corr;
        E_samp{ii}(jj)=tE;
        E0_samp{ii}(jj)=tE0;
    end
end

%%% statistics
g2_mean=cell(size(g2_bootstrap));
g2_sdev=cell(size(g2_bootstrap));
for ii=1:nparam
    g2_mean{ii}=cellfun(@(x) mean(x,4,'omitnan'),g2_bootstrap{ii},'UniformOutput',false);
    g2_sdev{ii}=cellfun(@(x) std(x,[],4,'omitnan'),g2_bootstrap{ii},'UniformOutput',false);
end

E_bootstrap_mean=cellfun(@(x) mean(x,'omitnan'),E_samp);
E_bootstrap_sdev=cellfun(@(x) std(x,'omitnan'),E_samp);

E0_bootstrap_mean=cellfun(@(x) mean(x,'omitnan'),E0_samp);
E0_bootstrap_sdev=cellfun(@(x) std(x,'omitnan'),E0_samp);


% % DEBUG suppressed
% %%% diplay some output
% par_T'
% 
% E_par'
% [E_bootstrap_mean,E_bootstrap_sdev]'
% 
% E0_par'
% [E0_bootstrap_mean,E0_bootstrap_sdev]'


%% PLOT: g2 visualisation
[cc,clight,cdark]=palette(3);   % colors
mark_type={'o','^','d'};        % markers for each axis
mark_size=7;
line_wid=1.1;


hg2_3d=[];
for ii=1:nparam
    hg2_3d(ii)=figure('Name',sprintf('g2_3d_T%0.3g',par_T(ii)));
    
    %%% spin-permutations as subplots
    perm_ord=[2,3,1];
    for jj=1:3
        subplot(1,3,jj);
        
        % permute 1D-slice through each Cart axis
        temp_ord=[1,2,3];
        temp_g2_perm=g2{ii}{jj};    % temporary var to hold dimension permuted g2
        temp_g2_sdev_perm=g2_sdev{ii}{jj};
        temp_g2_fit_perm=g2_fit{ii}{jj};
        
        for kk=1:3
            hold on;
            temp_ord=temp_ord(perm_ord);
            temp_g2_perm=permute(temp_g2_perm,perm_ord);
            temp_g2_sdev_perm=permute(temp_g2_sdev_perm,perm_ord);
            temp_g2_fit_perm=permute(temp_g2_fit_perm,perm_ord);
            
            % data
            th=ploterr(dk_cent_vec,temp_g2_perm(:,idx_dk0,idx_dk0),...
                [],temp_g2_sdev_perm(:,idx_dk0,idx_dk0),...
                mark_type{kk},'hhxy',0);
            set(th(1),'color',cc(kk,:),'Marker',mark_type{kk},'LineWidth',line_wid,...
                'MarkerSize',mark_size,'MarkerFaceColor',clight(kk,:));
            set(th(2),'color',cc(kk,:),'LineWidth',line_wid);
            
            % fitted model
            th=plot(dk_fit_vec,temp_g2_fit_perm(:,idx_dk0_fit,idx_dk0_fit),...
                'Color',clight(kk,:),'LineWidth',line_wid);
            uistack(th,'bottom');
        end
        
        % annotate
        ax=gca;
        ylim0=ax.YLim;
        ylim([0,ylim0(2)]);
        
        xlabel('$\Delta k$');
        ylabel('$g^{(2)}$');
        box on;
    end
end


%% DEBUG OUTPUT
for ii=1:3
    fprintf('%0.3f\n',g2_sdev{1}{ii}(idx_dk0,idx_dk0,idx_dk0));
end


%% Correlation - data vis
% Measured correlation
figure('Name','E_raw');

hh=ploterr(par_T/10e-6,E_par,...
    [],E_bootstrap_sdev,...
    'o','hhxy',0);
set(hh(1),'Marker','o','MarkerSize',7,...
    'Color','b','LineWidth',1.2,...
    'MarkerFaceColor','w');
set(hh(2),'Color','b','LineWidth',1.2);  % Y-err

% xlim([0,1]);
ylim([-1,1]);
xlabel('$\theta_0/\pi$');
ylabel('$E$');


% Mode occupancy corrected correlation
figure('Name','E_corrected');

hh=ploterr(par_T/10e-6,E0_par,...
    [],E0_bootstrap_sdev,...
    'o','hhxy',0);
set(hh(1),'Marker','o','MarkerSize',7,...
    'Color','k','LineWidth',1.2,...
    'MarkerFaceColor','w');
set(hh(2),'Color','k','LineWidth',1.2);  % Y-err

% xlim([0,1]);
ylim([-1,1]);
xlabel('$\theta_0/\pi$');
ylabel('$\bar{E}$');

ylim([-1,1]);



% % DEBUG below suppressed
% %% simple vis of 3D g2
% %   NOTE: 3d g2 is completely averaged across Y-axis
% for ii=1:nparam
%     
%     tT=1e6*par_T(ii);
%     tfigname=sprintf('g2_T%0.3g',tT);
%     
%     figure('Name',tfigname);
%     hold on;
%     for jj=1:3      % 3-types of spin-spin type of mom-correlation
%         subplot(1,3,jj);
%         imagesc(mean(g2{ii}{jj},3));        % averaged across Y (3rd dim)
%         
%         title(sprintf('TH=%d, G2=%d',ii,jj));
%         colorbar;
%         
%         axis square;
%     end
% end


%% g2 spatial distribution
% % data to analyse
% k=k_par{1};
% 
% % config
% naz=200;
% nel=100;
% [az,el]=sphgrid(naz,nel);
% 
% b_pole=(abs(el)>asin(0.6));      % bool to bad region around poles
% 
% 
% % count
% NN=cellfun(@(x) haloZoneCount(x,az,el,0.03,[],'simple'),k,'UniformOutput',false);
% Nhalo{1}=cat(3,NN{:,1});
% Nhalo{2}=cat(3,NN{:,2});
% 
% Nhalo{1}(repmat(b_pole,[1,1,nshot]))=NaN;
% Nhalo{2}(repmat(b_pole,[1,1,nshot]))=NaN;
% 
% % analysis
% summary_disthalo_g2dist(Nhalo,az,el,1);
% 
% 

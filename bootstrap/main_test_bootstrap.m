%% Test bootstrapping g2 analysis
%   [ ] convergence with more samples
%   [ ] independence of subset size
%   [ ] statistical variation
%
% DKS
% 20180811


%% config for bootstrapping
n_bs_rep=1;        % no. times to repeat bootstrapping for statistical variation
n_frac_samp_vec=[0.05,0.01];     % fractional size of sample
n_subset_vec=[10,20,40,80,160,320,640];         % no. of subsets


%% evaluate all configurations
niter_frac_samp=numel(n_frac_samp_vec);
niter_subset=numel(n_subset_vec);

g2_bs_avg=cell(niter_frac_samp,1);
g2_bs_std=cell(niter_frac_samp,1);
for idx_frac=1:niter_frac_samp
    g2_bs_avg{idx_frac}=NaN(niter_subset,3);
    g2_bs_std{idx_frac}=NaN(niter_subset,3);
    for idx_samp=1:niter_subset
        n_frac_samp=n_frac_samp_vec(idx_frac);    % subset fractional size
        n_subset=n_subset_vec(idx_samp);          % no. of subsets
        
        tg2_bs_sdev=NaN(n_bs_rep,3);
        for idx_rep=1:n_bs_rep
            % call bootstrapping script
            run('run_bootstrap.m');
            
            % get representative variability
            tg2_bs_sdev(idx_rep,:)=cellfun(@(x) x(idx_dk0,idx_dk0,idx_dk0),g2_sdev{1});
        end
        % statistics and store
        g2_bs_avg{idx_frac}(idx_samp,:)=mean(tg2_bs_sdev,1,'omitnan');
        g2_bs_std{idx_frac}(idx_samp,:)=std(tg2_bs_sdev,[],1,'omitnan');
    end
end


%% VIS
% config
[c0,clight,cdark]=palette(niter_frac_samp);
mark_type={'o','^','d','s'};    % markers for each axis
mark_size=7;
line_wid=1.1;


hfig=figure('Name','bootstrap_g2');

% data
p=NaN(niter_frac_samp ,1);
for idx_frac=1:niter_frac_samp 
    hold on;
    
    % just plot one of the g2 configurations
    idx_plot_g2=3;
    
    % variability
    tp=ploterr(n_subset_vec,g2_bs_avg{idx_frac}(:,idx_plot_g2),...
        [],g2_bs_std{idx_frac}(:,idx_plot_g2),...
        mark_type{mod(idx_frac,numel(mark_type))+1},'hhxy',0);
    set(tp(1),'color',c0(idx_frac,:),'Marker',mark_type{idx_frac},'LineWidth',line_wid,...
        'MarkerSize',mark_size,'MarkerFaceColor',clight(idx_frac,:),...
        'DisplayName',sprintf('%0.3g',n_frac_samp_vec(idx_frac)));
    set(tp(2),'color',c0(idx_frac,:),'LineWidth',line_wid);
    p(idx_frac)=tp(1); 
    
    % trend
    tline=plot(n_subset_vec,g2_bs_avg{idx_frac}(:,idx_plot_g2),...
        'Color',clight(idx_frac,:),'Marker','none','LineStyle','-',...
        'LineWidth',line_wid,...
        'DisplayName',sprintf('%0.3g',n_frac_samp_vec(idx_frac)));
    uistack(tline,'bottom');
end
% annotation
lgd=legend(p);
lgd.Title.String='Subset fraction';
xlabel('Number of samples');
ylabel('Error from bootstrapping');
box on;
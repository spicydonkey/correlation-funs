function [g2,G,N,h]=g2_ang_disthalo(k,dth_ed,bplot)
% G2_ANG_DISTHALO evaluates g2 correlation function in angular coordinates
% for 2-state resolved k-vectors
%
% INPUT
%   k: #shot x 2 cell array of zxy-vectors
%   dth_ed: vector of rel angle histogram bin edges
%   bplot: boolean flag for graphics
%
% OUTPUT
%   g2: 1x3 cell array of g2 with corr-types: (1,1), (2,2) and (1,2)
%   G: 1x3 cell array correlated # histogram (normalised)
%   N: 1x3 cell array uncorrelated # histogram (normalised)
%   h: output figure handle
%
% DKS
% 2018-10-27

%% g2
dth=edge2cent(dth_ed);
dth_deg=rad2deg(dth);

% preallocate
g2=cell(1,3);
G=cell(1,3);
N=cell(1,3);

% 11 / 22
for ii=1:2
    [g2{ii},G{ii},N{ii}]=g2_ang(k(:,ii),dth_ed);
end
% 12
[g2{3},G{3},N{3}]=g2x_ang(k,dth_ed);

%% vis
% configs
f_units='normalized';
% f_pos=[0.2,0.2,0.2,0.3];
f_pos=[0.2,0.2,0.4,0.6];
f_ren='painters';

[c,cl,cd]=palette(3);
line_sty={'-','--',':'};
mark_typ={'o','d','^'};
str_ss={'$(1,1)$','$(2,2)$','$(1,2)$'};
mark_siz=7;
line_wid=1.5;
fontsize=12;
ax_lwidth=1.2;

% figure
if bplot
    pline=NaN(3,1);
    h=figure('Units',f_units,'Position',f_pos,'Renderer',f_ren);
    hold on;
    for ii=1:3
        pline(ii)=plot(dth,g2{ii},'Color',c(ii,:),'LineStyle',line_sty{ii},...
            'LineWidth',line_wid,'MarkerSize',mark_siz,...
            'Marker',mark_typ{ii},'MarkerFaceColor',cl(ii,:),...
            'DisplayName',str_ss{ii});
    end
    % annotate
    box on;
    ax=gca;
    ax.LineWidth=ax_lwidth;
    ax.FontSize=fontsize;
    ax.LineWidth=1.2;
    xlabel('$\theta$ (rad)');
    ylabel('$g^{(2)}$');
    legend(pline);
end

end
function [ f ] = func_plot_slope_variances(struct_stat_beta, table_sign_beta, struct_stat_rsn_beta, chanlocs, vtype, frange, flag_sign)

% Function to show and contrast spectral slope (Beta) between the HC and SZ
% groups. Upper panels show channel-wise topology, lower panel shows
% statistics on the level of resting-state networks. Showing channel-wise
% differences is optional (flag_sign). Use vtype='mean' or vtype='var' for
% mean or variance over time, respectively. Use frange='lo' or frange='hi'
% for low- or high-frequency ranges, respectively.
%
% Note that this function uses the topoplot() function of EEGLAB (Delorme
% & Makeig, 2004) and therefore EEGLAB (or topoplot) needs to be added to
% the path.
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 02/07/2025


ind_beta = (contains(table_sign_beta.meas,[vtype '_Beta_' frange])==1); % select sata
p_beta = table_sign_beta(ind_beta,:);
f = figure('color','w','units','normalized','outerposition',[0.05 0.05 0.9 0.8]);

%% obtain channel-wise results
if strcmp(frange, 'hi')
    % ensure that data is sorted in original channel order
    table_beta = sortrows(struct_stat_beta.Beta_hi.var,'chID','ascend');
    cblabel = '\sigma(\beta_{\it{hi}})'; % label for colorbar
    str_title = '[20-45 Hz]'; % indicating frequency range in title
elseif strcmp(frange, 'lo')
    % ensure that data is sorted in original channel order
    table_beta = sortrows(struct_stat_beta.Beta_lo.var,'chID','ascend');
    cblabel = '\sigma(\beta_{\it{lo}})'; % label for colorbar
    str_title = '[1-4 Hz]'; % indicating frequency range in title
else
    error('Set frange to hi or lo.')
end

% get group expected values
v_hc = table_beta.E_hc;
v_sz = table_beta.E_sz;

% obtain global minimum and maximum
min_glob = min([v_hc; v_sz]);
max_glob = max([v_hc; v_sz]);

% denote significant channel-wise differences
ch_sign = p_beta.chID;

%% plotting

% HC topoplot
subplot(2,2,1)
if flag_sign % denote channel-wise between-group differences with asterisk
    topoplot(v_hc, chanlocs, 'emarker', {'.','k',20,2}, 'emarker2', {ch_sign,'*','k',20,3}, 'whitebk', 'on');
else 
    topoplot(v_hc, chanlocs, 'emarker', {'.','k',20,2}, 'whitebk', 'on');
end

colormap parula
clim([min_glob max_glob])
cb1 = colorbar;
cb1.Label.String = cblabel;
cb1.Label.Rotation = 270;
cb1.FontSize = 16;
cb1.Label.FontSize = 28;
cb1.Label.VerticalAlignment = 'bottom';
title(['Healthy, ' str_title],'FontSize',28,'visible','on')

% SZ topoplot
subplot(2,2,2)
if flag_sign % denote channel-wise between-group differences with asterisk
    topoplot(v_sz, chanlocs, 'emarker', {'.','k',20,2}, 'emarker2', {ch_sign,'*','k',20,3}, 'whitebk', 'on');
else
    topoplot(v_sz, chanlocs, 'emarker', {'.','k',20,2}, 'whitebk', 'on');
end

clim([min_glob max_glob])
colormap parula
cb2 = colorbar;
cb2.Label.String = cblabel;
cb2.Label.Rotation = 270;
cb2.FontSize = 16;
cb2.Label.FontSize = 28;
cb2.Label.VerticalAlignment = 'bottom';
title(['Schizophrenia, ' str_title],'FontSize',28,'visible','on')

%% RSN-wise results
if strcmp(frange,'hi')
    % ensure that data is sorted in the original order of RSNs
    rsn_beta = sortrows(struct_stat_rsn_beta.Beta_hi.var,'rsnID','ascend');
    str_ylabel = '\sigma(\beta_{\it{hi}})';
    str_title = '\beta_{hi}';
elseif strcmp(frange,'lo')
    % ensure that data is sorted in the original order of RSNs
    rsn_beta = sortrows(struct_stat_rsn_beta.Beta_lo.var,'rsnID','ascend');
    str_ylabel = '\sigma(\beta_{\it{lo}})';
    str_title = '\beta_{lo}';
end

rsn_list = rsn_beta.rsn;
nrsn = length(rsn_list);

% create box plots
subplot(2,2,[3,4])
for rsn = 1:nrsn
    v_hc = rsn_beta.v_hc{rsn};
    N_hc = length(v_hc);
    v_sz = rsn_beta.v_sz{rsn};
    N_sz = length(v_sz);

    % plot values
    plot((rsn-0.2)*ones(N_hc,1),v_hc,'o','color',[0,0.45,1],'MarkerSize',12)
    hold on
    plot((rsn+0.2)*ones(N_sz,1),v_sz,'o','color',[1,0.45,0],'MarkerSize',12)

    % draw boxes
    pm_hc = fill([rsn-0.05 rsn-0.35 rsn-0.35 rsn-0.05],[quantile(v_hc,0.25), quantile(v_hc,0.25), quantile(v_hc,0.75) quantile(v_hc,0.75)],[0 0.45 1],'LineWidth',1.5,'FaceAlpha',0.4);
    pm_sz = fill([rsn+0.05 rsn+0.35 rsn+0.35 rsn+0.05],[quantile(v_sz,0.25), quantile(v_sz,0.25), quantile(v_sz,0.75) quantile(v_sz,0.75)],[1 0.45 0],'LineWidth',1.5,'FaceAlpha',0.4);

    % draw percentiles (5, 25, 75, 95)
    plot([rsn-0.2 rsn-0.2],[quantile(v_hc,0.75) quantile(v_hc,0.95)],'k-','LineWidth',1.5)
    plot([rsn-0.3 rsn-0.1],[quantile(v_hc,0.95) quantile(v_hc,0.95)],'k-','LineWidth',1.5)
    plot([rsn-0.2 rsn-0.2],[quantile(v_hc,0.25) quantile(v_hc,0.05)],'k-','LineWidth',1.5)
    plot([rsn-0.3 rsn-0.1],[quantile(v_hc,0.05) quantile(v_hc,0.05)],'k-','LineWidth',1.5)

    plot([rsn+0.2 rsn+0.2],[quantile(v_sz,0.75) quantile(v_sz,0.95)],'k-','LineWidth',1.5)
    plot([rsn+0.1 rsn+0.3],[quantile(v_sz,0.95) quantile(v_sz,0.95)],'k-','LineWidth',1.5)
    plot([rsn+0.2 rsn+0.2],[quantile(v_sz,0.25) quantile(v_sz,0.05)],'k-','LineWidth',1.5)
    plot([rsn+0.1 rsn+0.3],[quantile(v_sz,0.05) quantile(v_sz,0.05)],'k-','LineWidth',1.5)

    % plot means and medians
    plot([rsn-0.35 rsn-0.05],[mean(v_hc) mean(v_hc)],'-','color',[0 0.35 0.8], 'LineWidth',3)
    plot([rsn-0.35 rsn-0.05],[median(v_hc) median(v_hc)],':','color',[0 0.35 0.8], 'LineWidth',3)
    plot([rsn+0.05 rsn+0.35],[mean(v_sz) mean(v_sz)],'-','color',[0.8 0.35 0], 'LineWidth',3)
    plot([rsn+0.05 rsn+0.35],[median(v_sz) median(v_sz)],':','color',[0.8 0.35 0], 'LineWidth',3)
end

ymax = max(get(gca,'ylim'));

% plot significant differences
for i = 1:nrsn
    p = rsn_beta.p(i);
    p_FDR = rsn_beta.p_FDR(i);
    if p_FDR < 0.05
        plot([i-0.3 i+0.3],[0.9*ymax 0.9*ymax],'k-','LineWidth',3)

        if (p < 0.05) && (p > 0.001)
            plot(i,0.95*ymax,'k*','MarkerSize',12,'LineWidth',2)
        elseif (p < 0.001) && (p > 0.0001)
            plot([i-0.025 i+0.025],[0.95*ymax 0.95*ymax],'k*','MarkerSize',12,'LineWidth',2)
        else
            plot([i-0.05 i i+0.05],[0.95*ymax 0.95*ymax 0.95*ymax],'k*','MarkerSize',12,'LineWidth',2)
        end
    elseif (p_FDR > 0.05) && (p < 0.05)
        plot([i-0.3 i+0.3],[0.9*ymax 0.9*ymax],':','color',[0.3,0.3,0.3],'LineWidth',1.5)
    end
end

% set panel parameters
xlim([0.5 nrsn+0.5])
set(gca,'LineWidth',1,'box','on','XTick',(1:nrsn),'XTickLabel',rsn_list,'FontSize',24);
legend([pm_hc,pm_sz],{'HC','SZ'},'FontSize',24,'location','northeast');
xlabel('Resting-state networks','FontSize',28)
ylabel(str_ylabel, 'FontSize', 28, 'visible', 'on')
title(['Temporal fluctuation of ' str_title ' in resting-state networks'],'FontSize',28,'visible','on')



end
function [ f ] = func_plot_slope_diff(table_hc, table_sz, chanlocs)

% Function to show the topological distribution of spectral indices for the
% low- and high-frequency regimes, as well as their difference (high-low).
%
% Note that this function uses the topoplot() function of EEGLAB (Delorme
% & Makeig, 2004) and therefore EEGLAB (or topoplot) needs to be added to
% the path.
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 02/07/2025

%% get data and constants

% ensure that data is sorted in original channel order
table_hc = sortrows(table_hc,'chID','ascend');
table_sz = sortrows(table_sz,'chID','ascend');

% find locations with significant differences
ch_sign_hc = find(table_hc.h_FDR);
ch_sign_sz = find(table_sz.h_FDR);

% find global minima and maxima for setting the color scale
min_glob = min([table_hc.E_lo; table_hc.E_hi; table_sz.E_lo; table_sz.E_hi]);
max_glob = max([table_hc.E_lo; table_hc.E_hi; table_sz.E_lo; table_sz.E_hi]);
min_diff = min([table_hc.E_hi - table_hc.E_lo; table_sz.E_hi - table_sz.E_lo]);
max_diff = max([table_hc.E_hi - table_hc.E_lo; table_sz.E_hi - table_sz.E_lo]);

%% plotting
f = figure('color','w','units','normalized','outerposition',[0 0 1 0.9]);

%% HC topoplots

% low-frequency regime (1-4 Hz)
subplot(2,3,1)
mu_lo = table_hc.E_lo;
topoplot(mu_lo, chanlocs, 'emarker', {'.','k',20,2}, 'whitebk', 'on');
clim([min_glob max_glob])
colormap parula
cb1 = colorbar;
cb1.Label.String = '\beta_{\it{lo}}';
cb1.Label.Rotation = 270;
cb1.FontSize = 16;
cb1.Label.FontSize = 28;
cb1.Label.VerticalAlignment = 'bottom';
title('\mu(\beta_{\it{lo}})','FontSize',28)
ylabel('Healthy','FontSize',28,'visible','on')

% high-frequency regime (20-45 Hz)
subplot(2,3,2)
mu_hi = table_hc.E_hi;
topoplot(mu_hi, chanlocs, 'emarker', {'.','k',20,2}, 'whitebk', 'on');
clim([min_glob max_glob])
colormap parula
cb2 = colorbar;
cb2.Label.String = '\beta_{\it{hi}}';
cb2.Label.Rotation = 270;
cb2.FontSize = 16;
cb2.Label.FontSize = 28;
cb2.Label.VerticalAlignment = 'bottom';
title('\mu(\beta_{\it{hi}})','FontSize',28)

% difference between low- and high-frequency regimes
subplot(2,3,3)
mu_diff = table_hc.E_hi - table_hc.E_lo;
topoplot(mu_diff, chanlocs, 'emarker',{'.','k',20,2},'emarker2',{ch_sign_hc,'*','k',16,3},'whitebk','on');
clim([min_diff max_diff])
colormap parula
cb3 = colorbar;
cb3.Label.String = '\Delta\beta';
cb3.Label.Rotation = 270;
cb3.FontSize = 16;
cb3.Label.FontSize = 28;
cb3.Label.VerticalAlignment = 'bottom';
title('\mu(\beta_{\it{hi}}-\beta_{\it{lo}})','FontSize',28)

%% SZ topoplots

% low-frequency regime (1-4 Hz)
subplot(2,3,4)
mu_lo = table_sz.E_lo;
topoplot(mu_lo, chanlocs, 'emarker', {'.','k',20,2}, 'whitebk', 'on');
clim([min_glob max_glob])
colormap parula
cb4 = colorbar;
cb4.Label.String = '\beta_{\it{lo}}';
cb4.Label.Rotation = 270;
cb4.FontSize = 16;
cb4.Label.FontSize = 28;
cb4.Label.VerticalAlignment = 'bottom';
title('\mu(\beta_{\it{lo}})','FontSize',28)
ylabel('Schizophrenia','FontSize',28,'visible','on')

% high-frequency regime (20-45 Hz)
subplot(2,3,5)
mu_hi = table_sz.E_hi;
topoplot(mu_hi, chanlocs, 'emarker', {'.','k',20,2}, 'whitebk', 'on');
clim([min_glob max_glob])
colormap parula
cb5 = colorbar;
cb5.Label.String = '\beta_{\it{hi}}';
cb5.Label.Rotation = 270;
cb5.FontSize = 16;
cb5.Label.FontSize = 28;
cb5.Label.VerticalAlignment = 'bottom';
title('\mu(\beta_{\it{hi}})','FontSize',28)

% difference between low- and high-frequency regimes
subplot(2,3,6)
mu_diff = table_sz.E_hi - table_sz.E_lo;
topoplot(mu_diff, chanlocs, 'emarker',{'.','k',20,2},'emarker2',{ch_sign_sz,'*','k',16,3},'whitebk','on');
clim([min_diff max_diff])
colormap parula
cb6 = colorbar;
cb6.Label.String = '\Delta\beta';
cb6.Label.Rotation = 270;
cb6.FontSize = 16;
cb6.Label.FontSize = 28;
cb6.Label.VerticalAlignment = 'bottom';
title('\mu(\beta_{\it{hi}}-\beta_{\it{lo}})','FontSize',28)

end
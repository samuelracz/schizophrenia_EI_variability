
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script first loads, structures and aggregates the outputs of the
% IRASA analysis, then performs the following statistical evaluations:
% - spectrum bimodality (independently for HC and SZ)
% - spectral slope mean over time
% - spectral slope variance over time
% - illustration of global EEG power spectra in HC and SZ
% - correlation analysis between findings and PANSS scores in SZ
% - demographics comparison between HC and SZ
%
% note that that is stored such that the first 30 subjects (sch_002 - 
% sch_050) are SZ patients, the remaining 31 subjects (sch_101 - sch_333)
% are healthy controls.
% 
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 02/06/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('functions'))
addpath('miscellaneous')

fnames_ec = dir('results_IRASA_ec/*.mat');
load('ws_55ch_montage.mat')

%% analysis constants

% map 55-channel EEG to resting-state networks (RSN)
Mr = func_map_to_RSN(M);

chlab = M.lab; %....................................channel labels
rsnlab = Mr.labels; %...............................RSN labels

ns = length(fnames_ec); %...total number of subjects
nch = length(chlab); %......total number of channels
nw = 12; %..................total number of sliding windows
n_sz = 30; %................number of schizophrenia patients
n_hc = 31; %................number of healthy controls
FDR_alpha = 0.05; %.........alpha for multiple comparisons adjustment

meanvars = {'mean','var'}; %........................analysis types
irasa_vars = {'Beta_lo','Beta_hi','Beta_bb'}; %.....frequency ranges
nrsn = length(rsnlab); %............................# of RSNs

%% preallocations

% structures for storing structured analysis outputs (ch and RSN level)
struct_hc = struct();
struct_sz = struct();
struct_hc_rsn = struct();
struct_sz_rsn = struct();

% structures for storing statistical testing outputs
struct_stat_beta = struct();
struct_stat_rsn_beta = struct();

for iv = 1:length(irasa_vars)
    struct_hc.(irasa_vars{iv}) = struct();
    struct_sz.(irasa_vars{iv}) = struct();
    for mv = 1:length(meanvars)
        % channel-level
        struct_hc.(irasa_vars{iv}).(meanvars{mv}) = zeros(n_hc, nch);
        struct_sz.(irasa_vars{iv}).(meanvars{mv}) = zeros(n_sz, nch);
        struct_stat_beta.(irasa_vars{iv}).(meanvars{mv}) = [];
        % RSN-level
        struct_hc_rsn.(irasa_vars{iv}).(meanvars{mv}) = zeros(n_hc, nrsn);
        struct_sz_rsn.(irasa_vars{iv}).(meanvars{mv}) = zeros(n_sz, nrsn);
        struct_stat_rsn_beta.(irasa_vars{iv}).(meanvars{mv}) = [];
    end
end

% for global spectra illustration
struct_hc.spec = struct('subjID',[],'plaw',[]);
struct_hc_rsn.spec = struct('subjID',[],'plaw',[]);
struct_sz.spec = struct('subjID',[],'plaw',[]);
struct_sz_rsn.spec = struct('subjID',[],'plaw',[]);

%% load and structure data

for subj = 1:ns % iterate over subjects
    load(['results_IRASA_ec/' fnames_ec(subj).name]);
    for iv = 1:length(irasa_vars) % iterate over fitting ranges
        v_tmp = zeros(nw, nch);
        for t = 1:nw % iterate over sliding windows
            v_tmp(t,:) = results_seg_ec(t).plaw.(irasa_vars{iv});
        end

        % aggregate over resting-state networks
        v_tmp_rsn = func_aggregate_RSNs_batch(v_tmp, Mr);

        if subj <= n_sz
            struct_sz.(irasa_vars{iv}).mean(subj,:) = mean(v_tmp,1);
            struct_sz.(irasa_vars{iv}).var(subj,:) = std(v_tmp,[],1);
            struct_sz_rsn.(irasa_vars{iv}).mean(subj,:) = mean(v_tmp_rsn,1);
            struct_sz_rsn.(irasa_vars{iv}).var(subj,:) = std(v_tmp_rsn,[],1);
        else
            struct_hc.(irasa_vars{iv}).mean(subj-n_sz,:) = mean(v_tmp,1);
            struct_hc.(irasa_vars{iv}).var(subj-n_sz,:) = std(v_tmp,[],1);
            struct_hc_rsn.(irasa_vars{iv}).mean(subj-n_sz,:) = mean(v_tmp_rsn,1);
            struct_hc_rsn.(irasa_vars{iv}).var(subj-n_sz,:) = std(v_tmp_rsn,[],1);
        end
    end

    % store data, compute average power spectrum for RSNs (illustration)
    if subj <= n_sz % first 30 subjects are SZ patients
        struct_sz.spec(subj).subjID = fnames_ec(subj).name;
        struct_sz_rsn.spec(subj).subjID = fnames_ec(subj).name;

        % compute average spectra (over sliding windows)
        [plaw, plaw_rsn] = func_average_RSN_results_irasa(results_seg_ec, chlab, Mr);
        struct_sz.spec(subj).plaw = plaw;
        struct_sz_rsn.spec(subj).plaw = plaw_rsn;
    else % remaining 31 subjects are HC individuals
        struct_hc.spec(subj-n_sz).subjID = fnames_ec(subj).name;
        struct_hc_rsn.spec(subj-n_sz).subjID = fnames_ec(subj).name;

        % compute average spectra (over sliding windows)
        [plaw, plaw_rsn] = func_average_RSN_results_irasa(results_seg_ec, chlab, Mr);
        struct_hc.spec(subj-n_sz).plaw = plaw;
        struct_hc_rsn.spec(subj-n_sz).plaw = plaw_rsn;
    end
end

% storing channel and RSN labels for convenience
struct_hc.chlab = chlab;
struct_sz.chlab = chlab;
struct_hc_rsn.rsnlab = rsnlab;
struct_sz_rsn.rsnlab = rsnlab;

%% check bimodality of the spectrum (channel-wise)

% mean of spectral slope over time (low- vs. high-range)
% channel-level
p_bdiff_mean_hc = func_pairwise_FDR(func_pairwise_bimodal_stat(struct_hc.Beta_lo.mean, struct_hc.Beta_hi.mean, chlab, 'Bdiff_mean'),0.05);
p_bdiff_mean_sz = func_pairwise_FDR(func_pairwise_bimodal_stat(struct_sz.Beta_lo.mean, struct_sz.Beta_hi.mean, chlab, 'Bdiff_mean'),0.05);

% RSN-level
p_bdiff_rsn_mean_hc = func_pairwise_FDR(func_pairwise_bimodal_stat(struct_hc_rsn.Beta_lo.mean, struct_hc_rsn.Beta_hi.mean, rsnlab, 'Bdiff_mean'),0.05);
p_bdiff_rsn_mean_sz = func_pairwise_FDR(func_pairwise_bimodal_stat(struct_sz_rsn.Beta_lo.mean, struct_sz_rsn.Beta_hi.mean, rsnlab, 'Bdiff_mean'),0.05);

% variance of spectral slope over time (low- vs. high-range - EXPLORATORY)
% channel-level
p_bdiff_var_hc = func_pairwise_FDR(func_pairwise_bimodal_stat(struct_hc.Beta_lo.var, struct_hc.Beta_hi.var, chlab, 'Bdiff_var'),0.05);
p_bdiff_var_sz = func_pairwise_FDR(func_pairwise_bimodal_stat(struct_sz.Beta_lo.var, struct_sz.Beta_hi.var, chlab, 'Bdiff_var'),0.05);

% RSN-level
p_bdiff_rsn_var_hc = func_pairwise_FDR(func_pairwise_bimodal_stat(struct_hc_rsn.Beta_lo.var, struct_hc_rsn.Beta_hi.var, rsnlab, 'Bdiff_var'),0.05);
p_bdiff_rsn_var_sz = func_pairwise_FDR(func_pairwise_bimodal_stat(struct_sz_rsn.Beta_lo.var, struct_sz_rsn.Beta_hi.var, rsnlab, 'Bdiff_var'),0.05);

%% run statistics - channel-wise (exploratory)

% preallocate
table_sign_beta = [];
table_sign_beta_unadjusted = []; % exploratory

for iv = 1:length(irasa_vars) % iterate over spectral indices
    for mv = 1:length(meanvars) % iterate over mean and variance

        % selecting variables
        v_hc = struct_hc.(irasa_vars{iv}).(meanvars{mv});
        v_sz = struct_sz.(irasa_vars{iv}).(meanvars{mv});
        vname = [meanvars{mv} '_' irasa_vars{iv}];
        
        % statistical testing
        p_corr = func_pairwise_FDR(func_pairwise_channel_stat(v_hc, v_sz, chlab, vname), FDR_alpha);

        % store outcomes
        struct_stat_beta.(irasa_vars{iv}).(meanvars{mv}) = p_corr;
    end
end

% collect significant differences
for iv = 1:length(irasa_vars)
    for mv = 1:length(meanvars)
        table_corr = struct_stat_beta.(irasa_vars{iv}).(meanvars{mv});
        tmp_sign = table_corr(table_corr.h==1,:);
        tmp_sign_FDR = table_corr(table_corr.h_FDR==1,:);

        % FDR-adjusted
        if isempty(tmp_sign_FDR)
            continue
        else
            if isempty(table_sign_beta)
                table_sign_beta = tmp_sign_FDR;
            else
                table_sign_beta = cat(1,table_sign_beta,tmp_sign_FDR);
            end
        end

        % unadjusted (exploratory)
        if isempty(tmp_sign)
            continue
        else
            if isempty(table_sign_beta_unadjusted)
                table_sign_beta_unadjusted = tmp_sign;
            else
                table_sign_beta_unadjusted = cat(1,table_sign_beta_unadjusted,tmp_sign);
            end
        end

        
    end
end



%% run statistics - RSN-wise

% preallocate
table_sign_rsn_beta = [];
table_sign_rsn_beta_unadjusted = []; % exploratory

for iv = 1:length(irasa_vars) % iterate over spectral indices
    for mv = 1:length(meanvars) % iterate over mean and variance

        % selecting variables
        v_hc = struct_hc_rsn.(irasa_vars{iv}).(meanvars{mv});
        v_sz = struct_sz_rsn.(irasa_vars{iv}).(meanvars{mv});
        vname = [meanvars{mv} '_' irasa_vars{iv}];

        % statistical testing
        p_corr = func_pairwise_FDR(func_pairwise_RSN_stat(v_hc, v_sz, rsnlab, vname), FDR_alpha);

        % store outcomes
        struct_stat_rsn_beta.(irasa_vars{iv}).(meanvars{mv}) = p_corr;
    end
end


% collect significant differences
for iv = 1:length(irasa_vars)
    for mv = 1:length(meanvars)
        table_corr = struct_stat_rsn_beta.(irasa_vars{iv}).(meanvars{mv});
        tmp_sign = table_corr(table_corr.h==1,:);
        tmp_sign_FDR = table_corr(table_corr.h_FDR==1,:);

        % FDR-adjusted
        if isempty(tmp_sign_FDR)
            continue
        else
            if isempty(table_sign_rsn_beta)
                table_sign_rsn_beta = tmp_sign_FDR;
            else
                table_sign_rsn_beta = cat(1,table_sign_rsn_beta,tmp_sign_FDR);
            end
        end

        % unadjusted (exploratory)
        if isempty(tmp_sign)
            continue
        else
            if isempty(table_sign_rsn_beta_unadjusted)
                table_sign_rsn_beta_unadjusted = tmp_sign;
            else
                table_sign_rsn_beta_unadjusted = cat(1,table_sign_rsn_beta_unadjusted,tmp_sign);
            end
        end
        
    end
end


%% saving structured data and between-group statistics
load('ws_55ch_chanlocs.mat')
save('miscellaneous/ws_results.mat', 'chanlocs', 'rsnlab', 'nch', 'nrsn',...
    'p_bdiff_mean_hc', 'p_bdiff_mean_sz', 'p_bdiff_var_hc', 'p_bdiff_var_sz',...
    'p_bdiff_rsn_mean_hc', 'p_bdiff_rsn_mean_sz', 'p_bdiff_rsn_var_hc', 'p_bdiff_rsn_var_sz',...
    'struct_hc', 'struct_sz', 'struct_stat_beta', ...
    'struct_hc_rsn', 'struct_sz_rsn', 'struct_stat_rsn_beta',...
    'table_sign_beta', 'table_sign_beta_unadjusted',...
    'table_sign_rsn_beta', 'table_sign_rsn_beta_unadjusted')

%% statistics for global spectra
[~, ~, spec_hc, spec_sz, struct_stat_spec] = func_plot_GA_spectra(struct_hc,struct_sz,0,0,0);

%% load and structure demographics and PANSS data

% HC demographics data
table_tmp_hc = readtable('demographics_table.xlsx','Sheet','HC');
table_dem_hc = table();
table_dem_hc.age = table_tmp_hc.age;
table_dem_hc.sex = table_tmp_hc.sex;
table_dem_hc.sexlab = table_tmp_hc.label;
table_dem_hc.edu_years = table_tmp_hc.edu_years;

% SZ demographics data
table_tmp_sz = readtable('demographics_table.xlsx','Sheet','SZ');
table_dem_sz = table();
table_dem_sz.age = table_tmp_sz.age;
table_dem_sz.sex = table_tmp_sz.sex;
table_dem_sz.sexlab = table_tmp_sz.label;
table_dem_sz.edu_years = table_tmp_sz.edu_years;
table_dem_sz.duration = table_tmp_sz.years;
table_dem_sz.CPZ = table_tmp_sz.CPZ;

% PANSS data
table_panss = readtable('demographics_table.xlsx','Sheet','SZ_correlation_analysis');

% define table of confounding variables
table_conf = table();
table_conf.age = table_panss.age;
table_conf.sex = table_panss.sex;
table_conf.edu = table_panss.edu_years;
table_conf.cpz = table_panss.CPZ;
table_conf.dur = table_panss.duration;
z_conf = table2array(table_conf);

%% demographics analyses

%% sex distributions
n1 = sum(table_dem_hc.sexlab);
N1 = length(table_dem_hc.sexlab);
n2 = sum(table_dem_sz.sexlab);
N2 = length(table_dem_sz.sexlab);

p0 = (n1+n2)/(N1+N2);
n10 = N1*p0;
n20 = N2*p0;

observed = [n1, N1-n1, n2, N2-n2];
expected = [n10, N1-n10, n20, N2-n20];

[~,p,stats] = chi2gof([1,2,3,4],'freq',observed,'expected',expected,'nparams',2);
disp(['sex: p=' num2str(p,'%.4f') ', chi-square=' num2str(stats.chi2stat,'%.4f')])

%% age
v_hc = table_dem_hc.age;
v_sz = table_dem_sz.age;

if lillietest(v_hc) || lillietest(v_sz)
    [p,~,stats] = ranksum(v_hc,v_sz);
    ttype = 'ranksum';
disp(['age: p=' num2str(p,'%.4f') ', ' ttype ', z=' num2str(stats.zval, '%.4f')])
else
    [~,p,~,stats] = ttest2(v_hc,v_sz);
    ttype = 'ttest2';
disp(['age: HC=' num2str(mean(v_hc),'%.2f') '+-' num2str(std(v_hc),'%.2f') ', SZ=' num2str(mean(v_sz),'%.2f') '+-' num2str(std(v_sz),'%.2f') ', p=' num2str(p,'%.4f') ', ' ttype ', t=' num2str(stats.tstat, '%.4f')])
end

%% years in education
v_hc = table_dem_hc.edu_years;
v_sz = table_dem_sz.edu_years;

if lillietest(v_hc) || lillietest(v_sz)
    [p,~,stats] = ranksum(v_hc,v_sz);
    ttype = 'ranksum';
    disp(['years in education: HC=' num2str(median(v_hc),'%.2f') ' [' num2str(min(v_hc)) '; ' num2str(max(v_hc)) '], SZ=' num2str(median(v_sz),'%.2f') ' [' num2str(min(v_sz)) '; ' num2str(max(v_sz)) '], p=' num2str(p,'%.4f') ', ' ttype ', z=' num2str(stats.zval, '%.4f')])
else
    [~,p,~,stats] = ttest2(v_hc,v_sz);
    ttype = 'ttest2';
disp(['years in education: p=' num2str(p,'%.4f') ', ' ttype ', t=' num2str(stats.tstat, '%.4f')])
end

disp(['illness duration: ' num2str(mean(table_dem_sz.duration),'%.2f') '+-' num2str(std(table_dem_sz.duration),'%.2f') ])
disp(['CPZ: ' num2str(mean(table_dem_sz.CPZ),'%.2f') '+-' num2str(std(table_dem_sz.CPZ), '%.2f') ])
disp(['PANSS SUM: ' num2str(mean(table_panss.SUM),'%.2f') '+-' num2str(std(table_panss.SUM), '%.2f') ])
disp(['PANSS GEN: ' num2str(mean(table_panss.GEN),'%.2f') '+-' num2str(std(table_panss.GEN), '%.2f') ])
disp(['PANSS NEG: ' num2str(mean(table_panss.NEG),'%.2f') '+-' num2str(std(table_panss.NEG), '%.2f') ])
disp(['PANSS POS: ' num2str(mean(table_panss.POS),'%.2f') '+-' num2str(std(table_panss.POS), '%.2f') ])

%% compute PANSS correlations - slope, channel-wise

corr_panss_beta_gen = func_corr_channel_PANSS(table_sign_beta, table_panss, table_conf, 'GEN');
corr_panss_beta_neg = func_corr_channel_PANSS(table_sign_beta, table_panss, table_conf, 'NEG');
corr_panss_beta_pos = func_corr_channel_PANSS(table_sign_beta, table_panss, table_conf, 'POS');
corr_panss_beta_sum = func_corr_channel_PANSS(table_sign_beta, table_panss, table_conf, 'SUM');

%% compute PANSS correlations - slope, RSN-wise

corr_panss_rsn_beta_gen = func_corr_RSN_PANSS(table_sign_rsn_beta, table_panss, table_conf, 'GEN');
corr_panss_rsn_beta_neg = func_corr_RSN_PANSS(table_sign_rsn_beta, table_panss, table_conf, 'NEG');
corr_panss_rsn_beta_pos = func_corr_RSN_PANSS(table_sign_rsn_beta, table_panss, table_conf, 'POS');
corr_panss_rsn_beta_sum = func_corr_RSN_PANSS(table_sign_rsn_beta, table_panss, table_conf, 'SUM');

%% saving results
load('ws_55ch_chanlocs.mat')
save('miscellaneous/ws_correlations.mat', 'chanlocs', 'rsnlab', 'nch', 'nrsn',...
    'table_panss', 'table_conf', 'z_conf',...
    'corr_panss_beta_gen', 'corr_panss_beta_neg', 'corr_panss_beta_pos', 'corr_panss_beta_sum',...
    'corr_panss_rsn_beta_gen', 'corr_panss_rsn_beta_neg', 'corr_panss_rsn_beta_pos', 'corr_panss_rsn_beta_sum')

%% generating statistics report table (optional)
% table_sign_write = func_generate_stat_report(table_sign);
% table_sign_hc_write = func_generate_stat_report_bimodal(p_bdiff_mean_hc(p_bdiff_mean_hc.h_FDR==1,:));
% table_sign_sz_write = func_generate_stat_report_bimodal(p_bdiff_mean_sz(p_bdiff_mean_sz.h_FDR==1,:));
% 
% writetable(table_sign_write,'table_stat_sign_v3.xlsx')
% writetable(table_sign_hc_write, 'table_stat_hc_sign_v3.xlsx')
% writetable(table_sign_sz_write, 'table_stat_sz_sign_v3.xlsx')



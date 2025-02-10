function [ fig1, fig2, output_hc, output_sz, struct_stat ] = func_plot_GA_spectra(struct_hc, struct_sz, ch, flag_fig, flag_slope)

% Function to visualize the mean and variance of power spectra (mixd, frac
% and osci) for both groups (HC and SZ). A specific channel can be selected
% (by channel index) or in case ch=0 the grand average spectrum over all
% channels is computed.
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 02/07/2025

%% preallocation and constants

% obtain spectra data
spec_hc = struct_hc.spec;
spec_sz = struct_sz.spec;

% select channel
if ch == 0 % global average over all channels
    chlab = 'global';
else % obtain channel label
    chlab = struct_hc.chlab{ch};
end

freq = spec_hc(1).plaw.freq; % frequency vector for spectra
N = length(freq); %............# of frequencies

n_hc = length(spec_hc); %......# of healthy controls
n_sz = length(spec_sz); %......# of SZ patients

% preallocations
hc_mixd = zeros(N,n_hc);
hc_frac = zeros(N,n_hc);
hc_osci = zeros(N,n_hc);

hc_mixdvar = zeros(N,n_hc);
hc_fracvar = zeros(N,n_hc);
hc_oscivar = zeros(N,n_hc);

sz_mixd = zeros(N,n_sz);
sz_frac = zeros(N,n_sz);
sz_osci = zeros(N,n_sz);
sz_mixdvar = zeros(N,n_sz);
sz_fracvar = zeros(N,n_sz);
sz_oscivar = zeros(N,n_sz);

%% analyses

% iterate over healthy controls
for n = 1:n_hc

    % get participant data
    tmp = spec_hc(n).plaw;
    if ch == 0 % global spectrum: average over channels
        hc_mixd(:,n) = mean(tmp.mixd,2);
        hc_frac(:,n) = mean(tmp.frac,2);
        hc_osci(:,n) = mean(tmp.osci,2);
        hc_mixdvar(:,n) = mean(tmp.mixdvar,2);
        hc_fracvar(:,n) = mean(tmp.fracvar,2);
        hc_oscivar(:,n) = mean(tmp.oscivar,2);
    else % specific channel: get channel spectra
        hc_mixd(:,n) = tmp.mixd(:,ch);
        hc_frac(:,n) = tmp.frac(:,ch);
        hc_osci(:,n) = tmp.osci(:,ch);
        hc_mixdvar(:,n) = tmp.mixdvar(:,ch);
        hc_fracvar(:,n) = tmp.fracvar(:,ch);
        hc_oscivar(:,n) = tmp.oscivar(:,ch);
    end
end

% iterate over SZ patients
for n = 1:n_sz

    % get participant data
    tmp = spec_sz(n).plaw;
    if ch == 0 % global spectrum: average over channels
        sz_mixd(:,n) = mean(tmp.mixd,2);
        sz_frac(:,n) = mean(tmp.frac,2);
        sz_osci(:,n) = mean(tmp.osci,2);
        sz_mixdvar(:,n) = mean(tmp.mixdvar,2);
        sz_fracvar(:,n) = mean(tmp.fracvar,2);
        sz_oscivar(:,n) = mean(tmp.oscivar,2);
    else % specific channel: get channel data
        sz_mixd(:,n) = tmp.mixd(:,ch);
        sz_frac(:,n) = tmp.frac(:,ch);
        sz_osci(:,n) = tmp.osci(:,ch);
        sz_mixdvar(:,n) = tmp.mixdvar(:,ch);
        sz_fracvar(:,n) = tmp.fracvar(:,ch);
        sz_oscivar(:,n) = tmp.oscivar(:,ch);
    end
end

%% storing output
output_hc = struct(...
    'freq', freq,...
    'mixd', hc_mixd,...
    'frac', hc_frac,...
    'osci', hc_osci,...
    'mixdvar', hc_mixdvar,...
    'fracvar', hc_fracvar,...
    'oscivar', hc_oscivar);

output_sz = struct(...
    'freq', freq,...
    'mixd', sz_mixd,...
    'frac', sz_frac,...
    'osci', sz_osci,...
    'mixdvar', sz_mixdvar,...
    'fracvar', sz_fracvar,...
    'oscivar', sz_oscivar);

%% resampling spectra for statistics (dimensionality reduction)
freq_res = 1:2:45;
freq_res_plot = 2:2:44;

hc_mixd_res = resample_spectrum(hc_mixd, freq, freq_res);
sz_mixd_res = resample_spectrum(sz_mixd, freq, freq_res);

hc_frac_res = resample_spectrum(hc_frac, freq, freq_res);
sz_frac_res = resample_spectrum(sz_frac, freq, freq_res);

hc_osci_res = resample_spectrum(hc_osci, freq, freq_res);
sz_osci_res = resample_spectrum(sz_osci, freq, freq_res);

hc_mixdvar_res = resample_spectrum(hc_mixdvar, freq, freq_res);
sz_mixdvar_res = resample_spectrum(sz_mixdvar, freq, freq_res);

hc_fracvar_res = resample_spectrum(hc_fracvar, freq, freq_res);
sz_fracvar_res = resample_spectrum(sz_fracvar, freq, freq_res);

hc_oscivar_res = resample_spectrum(hc_oscivar, freq, freq_res);
sz_oscivar_res = resample_spectrum(sz_oscivar, freq, freq_res);

%% computing statistics
p_mixd_spec = func_spectrum_stat(log(hc_mixd_res), log(sz_mixd_res), freq_res_plot, 0.05);
p_frac_spec = func_spectrum_stat(log(hc_frac_res), log(sz_frac_res), freq_res_plot, 0.05);
p_osci_spec = func_spectrum_stat(log(hc_osci_res), log(sz_osci_res), freq_res_plot, 0.05);

p_mixdvar_spec = func_spectrum_stat(log(hc_mixdvar_res), log(sz_mixdvar_res), freq_res_plot, 0.05);
p_fracvar_spec = func_spectrum_stat(log(hc_fracvar_res), log(sz_fracvar_res), freq_res_plot, 0.05);
p_oscivar_spec = func_spectrum_stat(log(hc_oscivar_res), log(sz_oscivar_res), freq_res_plot, 0.05);

%% within-group comparison of spectral power dynamics
[~, mh_mixd_hc, stat_mixd_hc] = func_spectrum_within_stat(hc_mixdvar_res,freq_res_plot,0.05);
[~, mh_osci_hc, stat_osci_hc] = func_spectrum_within_stat(hc_oscivar_res,freq_res_plot,0.05);

[~, mh_mixd_sz, stat_mixd_sz] = func_spectrum_within_stat(sz_mixdvar_res,freq_res_plot,0.05);
[~, mh_osci_sz, stat_osci_sz] = func_spectrum_within_stat(sz_oscivar_res,freq_res_plot,0.05);


%% storing output
struct_stat = struct(...
    'mixd', p_mixd_spec,...
    'frac', p_frac_spec,...
    'osci', p_osci_spec,...
    'mixdvar', p_mixdvar_spec,...
    'fracvar', p_fracvar_spec,...
    'oscivar', p_oscivar_spec,...
    'hc_mixd', stat_mixd_hc,...
    'hc_osci', stat_osci_hc,...
    'sz_mixd', stat_mixd_sz,...
    'sz_osci', stat_osci_sz);

%% plotting

if flag_fig
    fig1 = figure('color','w','units','normalized','outerposition',[0 0 1 0.9]);

    %% spectra mean over time 
    for i = 1:3
        switch i
            case 1 % mixed spectra
                v_hc = hc_mixd;
                v_sz = sz_mixd;
                stype = '\mu(BLP_{\it{mixd}})';
                xlab = 'Log(frequency)';
                ylab = 'Log(amplitude)';
                statvar = p_mixd_spec;
            case 2 % isolated fractal spectra
                v_hc = hc_frac;
                v_sz = sz_frac;
                stype = '\mu(BLP_{\it{frac}})';
                xlab = 'Log(frequency)';
                ylab = 'Log(amplitude)';
                statvar = p_frac_spec;
            case 3 % isolated oscillatory spectra
                v_hc = hc_osci;
                v_sz = sz_osci;
                stype = '\mu(BLP_{\it{osci}})';
                xlab = 'Frequency';
                ylab = 'Amplitude';
                statvar = p_osci_spec;
        end

        if i==1 || i==2
            mhc = mean(v_hc,2); % HC mean
            mhcp = mean(v_hc,2) + std(v_hc,[],2)./sqrt(n_hc); % HC mean + SEM
            mhcn = mean(v_hc,2) - std(v_hc,[],2)./sqrt(n_hc); % HC mean - SEM

            msz = mean(v_sz,2); % SZ mean
            mszp = mean(v_sz,2) + std(v_sz,[],2)./sqrt(n_sz); % SZ mean + SEM
            mszn = mean(v_sz,2) - std(v_sz,[],2)./sqrt(n_sz); % SZ mean - SEM
            
            % plotting
            subplot(2,3,i)
            phc = loglog(freq, mhc, 'b', 'LineWidth', 3);
            hold on
            psz = loglog(freq, msz, 'r', 'LineWidth', 3);
            xlim([1 45])
            
            % slopes
            if flag_slope
                % frequency range vectors
                f_lo = freq(freq>=1 & freq<=4);
                f_hi = freq(freq>=20 & freq<=45);

                % spectral powers
                s_lo_hc = mhc(freq>=1 & freq<=4);
                s_hi_hc = mhc(freq>=20 & freq<=45);
                s_lo_sz = msz(freq>=1 & freq<=4);
                s_hi_sz = msz(freq>=20 & freq<=45);

                % linear fit
                p_lo_hc = polyfit(log(f_lo),log(s_lo_hc),1);
                sl_lo_hc = exp(polyval(p_lo_hc,log(f_lo)));
                p_hi_hc = polyfit(log(f_hi),log(s_hi_hc),1);
                sl_hi_hc = exp(polyval(p_hi_hc,log(f_hi)));

                p_lo_sz = polyfit(log(f_lo),log(s_lo_sz),1);
                sl_lo_sz = exp(polyval(p_lo_sz,log(f_lo)));
                p_hi_sz = polyfit(log(f_hi),log(s_hi_sz),1);
                sl_hi_sz = exp(polyval(p_hi_sz,log(f_hi)));

                % adding slopes to the plot
                loglog(f_lo, sl_lo_hc,'--','color',[0.2,0.2,1],'LineWidth',4)
                loglog(f_hi, sl_hi_hc,'--','color',[0.2,0.2,1],'LineWidth',4)
                loglog(f_lo, sl_lo_sz,'--','color',[1,0.2,0.2],'LineWidth',4)
                loglog(f_hi, sl_hi_sz,'--','color',[1,0.2,0.2],'LineWidth',4)
            end
        else % i==3: plot isolated oscillatory spectra (on linear scale)
            mhc = mean(v_hc,2); % HC mean
            mhcp = mean(v_hc,2) + std(v_hc,[],2)./sqrt(n_hc); % HC mean + SEM
            mhcn = mean(v_hc,2) - std(v_hc,[],2)./sqrt(n_hc); % HC mean - SEM

            msz = mean(v_sz,2); % SZ mean
            mszp = mean(v_sz,2) + std(v_sz,[],2)./sqrt(n_sz); % SZ mean + SEM
            mszn = mean(v_sz,2) - std(v_sz,[],2)./sqrt(n_sz); % SZ mean - SEM

            % plot spectra
            subplot(2,3,i)
            phc = plot(freq, mhc, 'b', 'LineWidth', 3);
            hold on
            psz = plot(freq, msz, 'r', 'LineWidth', 3);
            xlim([1 45])
        end

        % error ranges
        fill([freq; flipud(freq)],[mhc; flipud(mhcp)],[0.5 0.5 1],'FaceAlpha',0.3,'EdgeAlpha',0.1)
        fill([freq; flipud(freq)],[mhc; flipud(mhcn)],[0.5 0.5 1],'FaceAlpha',0.3,'EdgeAlpha',0.1)
        fill([freq; flipud(freq)],[msz; flipud(mszp)],[1 0.5 0.5],'FaceAlpha',0.3,'EdgeAlpha',0.1)
        fill([freq; flipud(freq)],[msz; flipud(mszn)],[1 0.5 0.5],'FaceAlpha',0.3,'EdgeAlpha',0.1)

        % plotting statistics
        N = size(statvar,1);
        ymin = min(get(gca,'ylim'));
        ymax = max(get(gca,'ylim'));
        yrange = ymax - ymin;

        for n = 1:N
            if statvar.h_FDR(n) == 1
                if i == 1 || i == 2
                    plot([statvar.freq(n) statvar.freq(n)], ([ymax ymax]), 'k*', 'MarkerSize', 12, 'LineWidth', 2);
                    ylim([ymin, ymax+0.5*yrange])
                elseif i == 3
                    plot([statvar.freq(n) statvar.freq(n)], ([ymax ymax]), 'k*', 'MarkerSize', 12, 'LineWidth', 2);
                    ylim([ymin, ymax+0.1*yrange])
                end
            end
        end

        if i == 1 || i == 2
            ymin_tmp = min(get(gca,'ylim'));
            ymax_tmp = max(get(gca,'ylim'));
            plot([4 4],[ymin_tmp ymax_tmp],'k--','LineWidth',1.5)
            plot([20 20],[ymin_tmp ymax_tmp],'k--','LineWidth',1.5)
            ylim([ymin_tmp ymax_tmp])
        end

        % plot labels
        set(gca,'FontSize',18,'LineWidth',1.5)
        title([chlab ' ' stype],'FontSize',28)
        xlabel(xlab,'FontSize',24)
        ylabel(ylab,'FontSize',24)
        if i==1
            legend([phc,psz],{'HC','SZ'},'location','southwest','FontSize',24)
        elseif i==2
            legend([phc,psz],{'HC','SZ'},'location','southwest','FontSize',24)
        elseif i==3
            legend([phc,psz],{'HC','SZ'},'location','northeast','FontSize',24)
        end

    end

    %% spectra variance over time
    for i = 1:3
        switch i
            case 1 % mixed spectra
                v_hc = hc_mixdvar;
                v_sz = sz_mixdvar;
                stype = '\sigma(BLP_{\it{mixd}})';
                xlab = 'Log(frequency)';
                ylab = 'Log(amplitude)';
                statvar = p_mixdvar_spec;
            case 2 % isolated fractal spectra
                v_hc = hc_fracvar;
                v_sz = sz_fracvar;
                stype = '\sigma(BLP_{\it{frac}})';
                xlab = 'Log(frequency)';
                ylab = 'Log(amplitude)';
                statvar = p_fracvar_spec;
            case 3 % isolated oscillatory spectra
                v_hc = hc_oscivar;
                v_sz = sz_oscivar;
                stype = '\sigma(BLP_{\it{osci}})';
                xlab = 'Frequency';
                ylab = 'Log(amplitude)';
                statvar = p_oscivar_spec;
        end

        mhc = mean(v_hc,2); % HC mean
        mhcp = mean(v_hc,2) + std(v_hc,[],2)./sqrt(n_hc); % HC mean + SEM
        mhcn = mean(v_hc,2) - std(v_hc,[],2)./sqrt(n_hc); % HC mean - SEM

        msz = mean(v_sz,2); % SZ mean
        mszp = mean(v_sz,2) + std(v_sz,[],2)./sqrt(n_sz); % SZ mean + SEM
        mszn = mean(v_sz,2) - std(v_sz,[],2)./sqrt(n_sz); % SZ mean - SEM

        subplot(2,3,i+3)
        if i == 1 || i == 2
            phc = loglog(freq, mhc, 'b', 'LineWidth', 3);
            hold on
            psz = loglog(freq, msz, 'r', 'LineWidth', 3);
            xlim([1 45])
        elseif i == 3
            phc = plot(freq, mhc, 'b', 'LineWidth', 3);
            hold on
            psz = plot(freq, msz, 'r', 'LineWidth', 3);
            xlim([1 45])
        end

        % error ranges
        fill([freq; flipud(freq)],[mhc; flipud(mhcp)],[0.5 0.5 1],'FaceAlpha',0.3,'EdgeAlpha',0.1)
        fill([freq; flipud(freq)],[mhc; flipud(mhcn)],[0.5 0.5 1],'FaceAlpha',0.3,'EdgeAlpha',0.1)
        fill([freq; flipud(freq)],[msz; flipud(mszp)],[1 0.5 0.5],'FaceAlpha',0.3,'EdgeAlpha',0.1)
        fill([freq; flipud(freq)],[msz; flipud(mszn)],[1 0.5 0.5],'FaceAlpha',0.3,'EdgeAlpha',0.1)

        % plotting statistics
        N = size(statvar,1);
        ymin = min(get(gca,'ylim'));
        ymax = max(get(gca,'ylim'));
        yrange = ymax - ymin;

        for n = 1:N
            if statvar.h_FDR(n) == 1
                if i == 1 || i == 2
                    plot([statvar.freq(n) statvar.freq(n)], ([ymax ymax]), 'k*', 'MarkerSize', 12, 'LineWidth', 2);
                    ylim([ymin, ymax+0.1*yrange])
                elseif i == 3
                    plot([statvar.freq(n) statvar.freq(n)], ([ymax ymax]), 'k*', 'MarkerSize', 12, 'LineWidth', 2);
                    ylim([ymin, ymax+0.1*yrange])
                end
            end
        end

        % plot labels
        set(gca,'FontSize',18,'LineWidth',1.5)
        title([chlab ' ' stype],'FontSize',28)
        xlabel(xlab,'FontSize',24)
        ylabel(ylab,'FontSize',24)
        if i==1
            legend([phc,psz],{'HC','SZ'},'location','northwest','FontSize',24)
        elseif i==2
            legend([phc,psz],{'HC','SZ'},'location','southwest','FontSize',24)
        elseif i==3
            legend([phc,psz],{'HC','SZ'},'location','northeast','FontSize',24)
        end

    end
else
    fig1 = [];
end

%% plotting h matrices - comparing variance in power on a frequency-to-frequency basis for mixd and osci power
if flag_fig
    Nfig = length(freq_res_plot);
    fig2 = figure('color','w','units','normalized','outerposition',[0.1 0.1 0.8 0.8]);

    % HC, variance of mixd power
    subplot(221)
    imagesc(mh_mixd_hc);
    axis square
    set(gca, 'colormap', copper)
    set(gca,'FontSize',16)
    set(gca,'XTick',2:2:Nfig,'XTickLabel',freq_res_plot(2:2:end))
    set(gca,'YTick',1:2:Nfig-1,'YTickLabel',freq_res_plot(1:2:end-1))
    title('Healthy, \sigma(BLP_{\it{mixd}})','FontSize',24)
    xlabel('Frequency (Hz)','FontSize',20)
    ylabel('Frequency (Hz)','FontSize',20)
    clim([0 1])
    cb1 = colorbar;
    cb1.Label.String = 'significance (\it{h})';
    cb1.Label.Rotation = 270;
    cb1.Label.FontSize = 20;
    cb1.Label.VerticalAlignment = 'bottom';
    set(cb1,'XTick',[0 1])

    % SZ, variance of mixd power
    subplot(222)
    imagesc(mh_mixd_sz);
    axis square
    set(gca, 'colormap', copper)
    set(gca,'FontSize',16)
    set(gca,'XTick',2:2:Nfig,'XTickLabel',freq_res_plot(2:2:end))
    set(gca,'YTick',1:2:Nfig-1,'YTickLabel',freq_res_plot(1:2:end-1))
    title('Schizophrenia, \sigma(BLP_{\it{mixd}})','FontSize',24)
    xlabel('Frequency (Hz)','FontSize',20)
    ylabel('Frequency (Hz)','FontSize',20)
    clim([0 1])
    cb2 = colorbar;
    cb2.Label.String = 'significance (\it{h})';
    cb2.Label.Rotation = 270;
    cb2.Label.FontSize = 20;
    cb2.Label.VerticalAlignment = 'bottom';
    set(cb2,'XTick',[0 1])

    % HC, variance of osci power
    subplot(223)
    imagesc(mh_osci_hc);
    axis square
    set(gca, 'colormap', sky)
    set(gca,'FontSize',16)
    set(gca,'XTick',2:2:Nfig,'XTickLabel',freq_res_plot(2:2:end))
    set(gca,'YTick',1:2:Nfig-1,'YTickLabel',freq_res_plot(1:2:end-1))
    title('Healthy, \sigma(BLP_{\it{osci}})','FontSize',24)
    xlabel('Frequency (Hz)','FontSize',20)
    ylabel('Frequency (Hz)','FontSize',20)
    clim([0 1])
    cb3 = colorbar;
    cb3.Label.String = 'significance (\it{h})';
    cb3.Label.Rotation = 270;
    cb3.Label.FontSize = 20;
    cb3.Label.VerticalAlignment = 'bottom';
    set(cb3,'XTick',[0 1])

    % SZ, variance of osci power
    subplot(224)
    imagesc(mh_osci_sz);
    axis square
    set(gca, 'colormap', sky)
    set(gca,'FontSize',16)
    set(gca,'XTick',2:2:Nfig,'XTickLabel',freq_res_plot(2:2:end))
    set(gca,'YTick',1:2:Nfig-1,'YTickLabel',freq_res_plot(1:2:end-1))
    title('Schizophrenia, \sigma(BLP_{\it{osci}})','FontSize',24)
    xlabel('Frequency (Hz)','FontSize',20)
    ylabel('Frequency (Hz)','FontSize',20)
    clim([0 1])
    cb4 = colorbar;
    cb4.Label.String = 'significance (\it{h})';
    cb4.Label.Rotation = 270;
    cb4.Label.FontSize = 20;
    cb4.Label.VerticalAlignment = 'bottom';
    set(cb4,'XTick',[0 1])

else
    fig2 = [];
end

end

%% helper functions
function [spec_res] = resample_spectrum(spec_in, freq_in, freq_res)
[~, Ns] = size(spec_in);
Nf = length(freq_res);
spec_res = zeros(Nf-1,Ns);

for f = 1:Nf-1
    spec_res(f,:) = mean(spec_in(freq_in>=freq_res(f) & freq_in<freq_res(f+1),:),1);
end
end

function [ table_out ] = func_spectrum_stat(data_hc, data_sz, freq, FDR_alpha)
N = length(freq);

table_out = table(...,
    zeros(N,1), zeros(N,1), zeros(N,1),...
    zeros(N,1), zeros(N,1), cell(N,1),...
    zeros(N,1), zeros(N,1),...
    'VariableNames', {'freq','E_hc','E_sz','p','h','ttype','p_FDR','h_FDR'});

% pairwise comparisons
for n = 1:N
    table_out.freq(n) = freq(n);
    v_hc = data_hc(n,:);
    v_sz = data_sz(n,:);
    if lillietest(v_hc) || lillietest(v_sz)
        [p,h] = ranksum(v_hc,v_sz);
        table_out.E_hc(n) = median(v_hc);
        table_out.E_sz(n) = median(v_sz);
        table_out.ttype{n} = 'ranksum';
    else
        [h,p] = ttest2(v_hc,v_sz);
        table_out.E_hc(n) = mean(v_hc);
        table_out.E_sz(n) = mean(v_sz);
        table_out.ttype{n} = 'ttest2';
    end
    table_out.p(n) = p;
    table_out.h(n) = h;
end

% multiple comparisons adjustment
table_out = sortrows(table_out,'p','ascend');

for n = 1:N
    table_out.p_FDR(n) = table_out.p(n)*N/n;
end

for n = 1:N
    if table_out.p_FDR(n) >= FDR_alpha
        break
    else
        table_out.h_FDR(n) = 1;
    end
end

% sorting to frequencies
table_out = sortrows(table_out,'freq','ascend');

end

function [mat_p, mat_h, table_out] = func_spectrum_within_stat(data, freq, FDR_alpha)
% within-group comparison of spectra

nf = length(freq);
struct_tmp = struct('f1',[],'f2',[],'E1',[],'E2',[],'p',[],'h',[],'p_FDR',[],'h_FDR',[],'ttype',[]);
ind_tmp = 1;
for f1 = 1:nf-1
    for f2 = f1+1:nf
        v1 = data(f1,:);
        v2 = data(f2,:);
        if lillietest(v1) || lillietest(v2)
            [p,h] = signrank(v1,v2);
            struct_tmp(ind_tmp).E1 = median(v1);
            struct_tmp(ind_tmp).E2 = median(v2);
            struct_tmp(ind_tmp).ttype = 'signrank';
        else
            [h,p] = ttest(v1,v2);
            struct_tmp(ind_tmp).E1 = mean(v1);
            struct_tmp(ind_tmp).E2 = mean(v2);
            struct_tmp(ind_tmp).ttype = 'ttest';
        end
        struct_tmp(ind_tmp).f1 = f1;
        struct_tmp(ind_tmp).f2 = f2;
        struct_tmp(ind_tmp).p = p;
        struct_tmp(ind_tmp).h = h;
        struct_tmp(ind_tmp).p_FDR = p;
        struct_tmp(ind_tmp).h_FDR = 0;
        ind_tmp = ind_tmp + 1;
    end
end

table_out = struct2table(struct_tmp);
% multiple comparisons adjustment
table_out = sortrows(table_out,'p','ascend');
N = size(table_out,1);

for n = 1:N
    table_out.p_FDR(n) = table_out.p(n)*N/n;
end

for n = 1:N
    if table_out.p_FDR(n) >= FDR_alpha
        break
    else
        table_out.h_FDR(n) = 1;
    end
end


% reconstructing matrices
mat_h = zeros(nf,nf);
mat_p = zeros(nf,nf);

for n = 1:N
    if table_out.h_FDR(n)
        f1 = table_out.f1(n);
        f2 = table_out.f2(n);
        mat_h(f1,f2) = 1;
        mat_h(f2,f1) = 1;
        mat_p(f1,f2) = table_out.p_FDR(n);
        mat_p(f2,f1) = table_out.p_FDR(n);
    end
end

end

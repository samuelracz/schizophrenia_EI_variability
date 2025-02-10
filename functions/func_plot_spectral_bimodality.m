function [fig1, output_hc, output_sz] = func_plot_spectral_bimodality(struct_hc, struct_sz, chlab1, chlab2, flag_fig)

% Function to illustrate spectral bimodality over two selected channels in
% HC and SZ groups. Channel selection is via channel labels.
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 02/07/2025

%% variables and constants

spec_hc = struct_hc.spec; %............HC spectra
spec_sz = struct_sz.spec; %............SZ spectra

chlab = struct_hc.chlab; %.............channel labels
ch1 = find(strcmp(chlab,chlab1)); %....channel 1 index
ch2 = find(strcmp(chlab,chlab2)); %....channel 2 index

freq = spec_hc(1).plaw.freq; %.........frequency vector
N = length(freq); %....................# of frequencies

n_hc = length(spec_hc); %..............# of healthy controls
n_sz = length(spec_sz); %..............# of szhizophrenia patients

%% obtain HC data

% preallocations
hc_mixd_ch1 = zeros(N,n_hc);
hc_frac_ch1 = zeros(N,n_hc);
hc_osci_ch1 = zeros(N,n_hc);
hc_mixd_ch2 = zeros(N,n_hc);
hc_frac_ch2 = zeros(N,n_hc);
hc_osci_ch2 = zeros(N,n_hc);

for n = 1:n_hc % iterate over HC participants
    tmp = spec_hc(n).plaw;
    hc_mixd_ch1(:,n) = tmp.mixd(:,ch1);
    hc_frac_ch1(:,n) = tmp.frac(:,ch1);
    hc_osci_ch1(:,n) = (log(tmp.mixd(:,ch1)) - log(tmp.frac(:,ch1)));
    hc_mixd_ch2(:,n) = tmp.mixd(:,ch2);
    hc_frac_ch2(:,n) = tmp.frac(:,ch2);
    hc_osci_ch2(:,n) = (log(tmp.mixd(:,ch2)) - log(tmp.frac(:,ch2)));
end

%% obtain SZ data

% preallocations
sz_mixd_ch1 = zeros(N,n_sz);
sz_frac_ch1 = zeros(N,n_sz);
sz_osci_ch1 = zeros(N,n_sz);
sz_mixd_ch2 = zeros(N,n_sz);
sz_frac_ch2 = zeros(N,n_sz);
sz_osci_ch2 = zeros(N,n_sz);

for n = 1:n_sz % iterate over SZ participants
    tmp = spec_sz(n).plaw;
    sz_mixd_ch1(:,n) = tmp.mixd(:,ch1);
    sz_frac_ch1(:,n) = tmp.frac(:,ch1);
    sz_osci_ch1(:,n) = (log(tmp.mixd(:,ch1)) - log(tmp.frac(:,ch1)));
    sz_mixd_ch2(:,n) = tmp.mixd(:,ch2);
    sz_frac_ch2(:,n) = tmp.frac(:,ch2);
    sz_osci_ch2(:,n) = (log(tmp.mixd(:,ch2)) - log(tmp.frac(:,ch2)));
end

%% store output
output_hc = struct(...
    'freq', freq,...
    'mixd_ch1', hc_mixd_ch1,...
    'frac_ch1', hc_frac_ch1,...
    'osci_ch1', hc_osci_ch1,...
    'mixd_ch2', hc_mixd_ch2,...
    'frac_ch2', hc_frac_ch2,...
    'osci_ch2', hc_osci_ch2);

output_sz = struct(...
    'freq', freq,...
    'mixd_ch1', sz_mixd_ch1,...
    'frac_ch1', sz_frac_ch1,...
    'osci_ch1', sz_osci_ch1,...
    'mixd_ch2', sz_mixd_ch2,...
    'frac_ch2', sz_frac_ch2,...
    'osci_ch2', sz_osci_ch2);



%% plotting 

xlims = [1 45];

if flag_fig
    fig1 = figure('color','w','units','normalized','outerposition',[0.2 0.05 0.6 0.9]);

    %% ploting spectra from HC

    for i = 1:2 % iterate over channels
        switch i
            case 1
                v_mixd = hc_mixd_ch1;
                v_frac = hc_frac_ch1;
                v_osci = hc_osci_ch1;
                xlab = 'Log(frequency)';
                ylab = 'Log(amplitude)';
            case 2
                v_mixd = hc_mixd_ch2;
                v_frac = hc_frac_ch2;
                v_osci = hc_osci_ch2;
                xlab = 'Log(frequency)';
                ylab = 'Log(amplitude)';
        end
        
        m_mixd = mean(v_mixd,2); % mean raw (mixd) power spectrum
        m_mixd_p = m_mixd + std(v_mixd,[],2)./sqrt(n_hc); % mean + SEM
        m_mixd_n = m_mixd - std(v_mixd,[],2)./sqrt(n_hc); % mean - SEM 

        m_frac = mean(v_frac,2); % mean fractal power spectrum
        % SEM range is not plotted for fractal power spectrum for clarity

        % comupte linear fits (slopes)
        f_lo = freq(freq>=1 & freq<=4);
        f_hi = freq(freq>=20 & freq<=45);
        s_lo = m_frac(freq>=1 & freq<=4);
        s_hi = m_frac(freq>=20 & freq<=45);

        p_lo_frac = polyfit(log(f_lo),log(s_lo),1);
        sl_lo_frac = exp(polyval(p_lo_frac,log(f_lo)));
        p_hi_frac = polyfit(log(f_hi),log(s_hi),1);
        sl_hi_frac = exp(polyval(p_hi_frac,log(f_hi)));

        % main axis for mixd and frac spectra
        ax1 = subplot(2,2,i,'Position',[(i-1)*0.5+0.08 0.58, 0.34, 0.36]);

        % error ranges
        fill([freq; flipud(freq)],[m_mixd_n; flipud(m_mixd_p)],[0.6 0.6 0.6],'FaceAlpha',0.3,'EdgeAlpha',0.1)
        xlim(xlims)
        set(ax1,'FontSize',18,'LineWidth',1.5,'XScale','log','YScale','log')
        hold on
        
        % plotting spectra
        pm = loglog(freq, m_mixd, 'color', [0.6, 0.6, 0.6], 'LineWidth', 3);
        pf = loglog(freq, m_frac, 'k', 'LineWidth', 3);

        % slopes
        loglog(f_lo, sl_lo_frac,'--','color',[0.2,0.2,1],'LineWidth',4)
        loglog(f_hi, sl_hi_frac,'--','color',[1,0.2,0.2],'LineWidth',4)

        % fit boundary
        ymin_tmp = min(get(gca,'ylim'));
        ymax_tmp = max(get(gca,'ylim'));
        plot([4 4],[ymin_tmp ymax_tmp],'k--','LineWidth',1.5)
        plot([20 20],[ymin_tmp ymax_tmp],'k--','LineWidth',1.5)
        ylim([ymin_tmp ymax_tmp])

        % shade boundary
        fill([4, 20, 20, 4], [ymin_tmp, ymin_tmp, ymax_tmp, ymax_tmp],[0.95,0.9,0.65],'FaceAlpha',0.2,'EdgeAlpha',0);

        % legend
        legend([pm,pf],{'mixed','fractal'},'location','northeast','FontSize',20)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NOTE: annotations are hard-coded to fit figures panels from F5
        % and P3, however ideal locations for showing \beta might change
        % for different channels

        % annotations       
        a_lo = annotation('textbox');
        a_lo.String = ['\beta_{lo}=' num2str(-p_lo_frac(1),'%.2f')];
        a_lo.FontSize = 24;
        a_lo.Color = 'b';
        if i == 1
            a_lo.Position = [(i-1)*0.5+0.1, 0.75, 0.1, 0.12];
        elseif i == 2
            a_lo.Position = [(i-1)*0.5+0.1, 0.7, 0.1, 0.12];
        end
        a_lo.LineStyle = 'none';

        a_hi = annotation('textbox');
        a_hi.String = ['\beta_{hi}=' num2str(-p_hi_frac(1),'%.2f')];
        a_hi.FontSize = 24;
        a_hi.Color = 'r';
        if i == 1
            a_hi.Position = [(i-1)*0.5+0.325, 0.62, 0.1, 0.12];
        elseif i == 2
            a_hi.Position = [(i-1)*0.5+0.325, 0.7, 0.1, 0.12];
        end
        a_hi.LineStyle = 'none';

        % plot labels
        if i == 1
            title([chlab1 ', Healthy'],'FontSize',28)
        elseif i == 2
            title([chlab2 ', Healthy'],'FontSize',28)
        end

        xlabel(xlab,'FontSize',24)
        ylabel(ylab,'FontSize',24)
        
        % inset plot for oscillatory spectrum
        ax2 = axes('Position',[(i-1)*0.5+0.1, 0.6, 0.1, 0.12]);

        m_osci = mean(v_osci,2); % mean oscillatory power
        m_osci_p = mean(v_osci,2) + std(v_osci,[],2)./sqrt(n_hc); % mean + SEM
        m_osci_n = mean(v_osci,2) - std(v_osci,[],2)./sqrt(n_hc); % mean - SEM

        % error ranges
        fill([freq; flipud(freq)],[m_osci_n; flipud(m_osci_p)],[0.6 0.6 0.6],'FaceAlpha',0.3,'EdgeAlpha',0.1)
        hold on

        % plot mean
        plot(freq, m_osci, 'k', 'LineWidth', 1);
        xlim(xlims)

        set(ax2,'FontSize',12)

    end

    %% plotting spectra from SZ

    for i = 1:2 % iterate over channels
        switch i
            case 1
                v_mixd = sz_mixd_ch1;
                v_frac = sz_frac_ch1;
                v_osci = sz_osci_ch1;
                xlab = 'Log(frequency)';
                ylab = 'Log(amplitude)';
            case 2
                v_mixd = sz_mixd_ch2;
                v_frac = sz_frac_ch2;
                v_osci = sz_osci_ch2;
                xlab = 'Log(frequency)';
                ylab = 'Log(amplitude)';
        end

        m_mixd = mean(v_mixd,2); % mean raw (mixd) power spectrum
        m_mixd_p = m_mixd + std(v_mixd,[],2)./sqrt(n_hc); % mean + SEM
        m_mixd_n = m_mixd - std(v_mixd,[],2)./sqrt(n_hc); % mean - SEM 

        m_frac = mean(v_frac,2); % mean fractal power spectrum
        % SEM range is not plotted for fractal power spectrum for clarity

        % slopes
        f_lo = freq(freq>=1 & freq<=4);
        f_hi = freq(freq>=20 & freq<=45);
        s_lo = m_frac(freq>=1 & freq<=4);
        s_hi = m_frac(freq>=20 & freq<=45);

        p_lo_frac = polyfit(log(f_lo),log(s_lo),1);
        sl_lo_frac = exp(polyval(p_lo_frac,log(f_lo)));
        p_hi_frac = polyfit(log(f_hi),log(s_hi),1);
        sl_hi_frac = exp(polyval(p_hi_frac,log(f_hi)));

        % main axis for mixd and frac spectra
        ax1 = subplot(2,2,i+2,'Position',[(i-1)*0.5+0.08 0.08, 0.34, 0.36]);

        % error ranges
        fill([freq; flipud(freq)],[m_mixd_n; flipud(m_mixd_p)],[0.3 0.3 0.3],'FaceAlpha',0.3,'EdgeAlpha',0.1)
        xlim(xlims)
        set(ax1,'FontSize',18,'LineWidth',1.5,'XScale','log','YScale','log')
        hold on
        
        % plotting spectra
        pm = loglog(freq, m_mixd, 'color', [0.6, 0.6, 0.6], 'LineWidth', 3);
        pf = loglog(freq, m_frac, 'k', 'LineWidth', 3);

        % slopes
        loglog(f_lo, sl_lo_frac,'--','color',[0.2,0.2,1],'LineWidth',4)
        loglog(f_hi, sl_hi_frac,'--','color',[1,0.2,0.2],'LineWidth',4)

        % fit boundary
        ymin_tmp = min(get(gca,'ylim'));
        ymax_tmp = max(get(gca,'ylim'));
        plot([4 4],[ymin_tmp ymax_tmp],'k--','LineWidth',1.5)
        plot([20 20],[ymax_tmp ymin_tmp],'k--','LineWidth',1.5)
        ylim([ymin_tmp ymax_tmp])

        % shade boundary
        fill([4, 20, 20, 4], [ymin_tmp, ymin_tmp, ymax_tmp, ymax_tmp],[0.95,0.9,0.65],'FaceAlpha',0.2,'EdgeAlpha',0);

        % legend
        legend([pm,pf],{'mixed','fractal'},'location','northeast','FontSize',20)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NOTE: annotations are hard-coded to fit figures panels from F5
        % and P3, however ideal locations for showing \beta might change
        % for different channels

        % annotations
        a_lo = annotation('textbox');
        a_lo.String = ['\beta_{lo}=' num2str(-p_lo_frac(1),'%.2f')];
        a_lo.FontSize = 24;
        a_lo.Color = 'b';
        if i == 1
            a_lo.Position = [(i-1)*0.5+0.1, 0.2, 0.1, 0.12];
        elseif i ==2
            a_lo.Position = [(i-1)*0.5+0.1, 0.225, 0.1, 0.12];
        end
        a_lo.LineStyle = 'none';

        a_hi = annotation('textbox');
        a_hi.String = ['\beta_{hi}=' num2str(-p_hi_frac(1),'%.2f')];
        a_hi.FontSize = 24;
        a_hi.Color = 'r';
        if i == 1
            a_hi.Position = [(i-1)*0.5+0.325, 0.15, 0.1, 0.12];
        elseif i == 2
            a_hi.Position = [(i-1)*0.5+0.325, 0.22, 0.1, 0.12];
        end
        a_hi.LineStyle = 'none';

        % plot labels
        if i == 1
            title([chlab1 ', Schizophrenia '],'FontSize',28)
        elseif i == 2
            title([chlab2 ', Schizophrenia '],'FontSize',28)
        end

        xlabel(xlab,'FontSize',24)
        ylabel(ylab,'FontSize',24)
        
        % inset plot for oscillatory spectrum
        ax2 = axes('Position',[(i-1)*0.5+0.1, 0.1, 0.1, 0.12]);

        m_osci = mean(v_osci,2); % mean oscillatory power
        m_osci_p = mean(v_osci,2) + std(v_osci,[],2)./sqrt(n_hc); % mean + SEM
        m_osci_n = mean(v_osci,2) - std(v_osci,[],2)./sqrt(n_hc); % mean - SEM

        % error ranges
        fill([freq; flipud(freq)],[m_osci_n; flipud(m_osci_p)],[0.6 0.6 0.6],'FaceAlpha',0.3,'EdgeAlpha',0.1)
        hold on

        % plot mean
        plot(freq, m_osci, 'k', 'LineWidth', 1);
        xlim(xlims)

        set(ax2,'FontSize',12)

    end
else
    fig1 = [];
end


end


function [ output_table ] = func_pairwise_channel_stat(data_hc, data_sz, chlab, meas)

% Function compare spectral indices obtained from HC and SZ groups
% (independent design) on the level of individual channels.
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 02/06/2025

%% preallocation and analysis constants

nch = length(chlab); % # of channels

output_table = table(cell(nch,1), cell(nch,1), zeros(nch,1),...
    zeros(nch,1), zeros(nch,1),...
    zeros(nch,1), zeros(nch,1),...
    cell(nch,1), zeros(nch,1),...
    zeros(nch,1), zeros(nch,1),...
    cell(nch,1), cell(nch,1),...
    'VariableNames', {'meas','ch', 'chID', 'E_hc','E_sz','p','h','ttype','tstat','p_FDR','h_FDR','v_hc','v_sz'});

%% statistical tests

for ch = 1:nch % iterate over channels

    % select variables
    v_hc = data_hc(:,ch);
    v_sz = data_sz(:,ch);

    % test for normality
    if lillietest(v_hc) || lillietest(v_sz)
        % nonparametric testing (Wilcoxon ranksum test)
        [p,h,s] = ranksum(v_hc, v_sz);
        E_hc = median(v_hc);
        E_sz = median(v_sz);
        ttype = 'ranksum';
        statvalue = s.zval;
    else
        % parametric testing (two sample t-test)
        [h,p,~,s] = ttest2(v_hc, v_sz);
        E_hc = mean(v_hc);
        E_sz = mean(v_sz);
        ttype = 'ttest2';
        statvalue = s.tstat;
    end

    % store outcomes
    output_table.meas{ch} = meas;
    output_table.ch{ch} = chlab{ch};
    output_table.chID(ch) = ch;
    output_table.E_hc(ch) = E_hc;
    output_table.E_sz(ch) = E_sz;
    output_table.ttype{ch} = ttype;
    output_table.tstat(ch) = statvalue;
    output_table.p(ch) = p;
    output_table.h(ch) = h;
    output_table.p_FDR(ch) = p;
    output_table.h_FDR(ch) = 0;
    output_table.v_hc{ch} = v_hc;
    output_table.v_sz{ch} = v_sz;
end

% sorting outcomes in descending order or significance (1-p)
output_table = sortrows(output_table, 'p', 'ascend');

end
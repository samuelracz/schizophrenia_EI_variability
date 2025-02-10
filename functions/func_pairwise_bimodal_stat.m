function [ output_table ] = func_pairwise_bimodal_stat(data_lo, data_hi, chlab, meas)

% Function compare spectral indices obtained from low- and high-frequency
% regimes (paired comparison design)
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
    'VariableNames', {'meas','ch', 'chID', 'E_lo','E_hi','p','h','ttype','tstat','p_FDR','h_FDR','v_lo','v_hi'});

%% statistical tests

for ch = 1:nch % iterate over channels

    % select variables
    v_lo = data_lo(:,ch);
    v_hi = data_hi(:,ch);

    % test for normality
    if lillietest(v_lo) || lillietest(v_hi)
        % nonparametric comparison - Wilcoxon signed rank test
        [p,h,s] = signrank(v_lo, v_hi);
        E_lo = median(v_lo);
        E_hi = median(v_hi);
        ttype = 'signrank';
        statvalue = s.zval;
    else
        % parametric comparison - paired t-test
        [h,p,~,s] = ttest(v_lo, v_hi);
        E_lo = mean(v_lo);
        E_hi = mean(v_hi);
        ttype = 'ttest';
        statvalue = s.tstat;
    end

    % store outcomes
    output_table.meas{ch} = meas;
    output_table.ch{ch} = chlab{ch};
    output_table.chID(ch) = ch;
    output_table.E_lo(ch) = E_lo;
    output_table.E_hi(ch) = E_hi;
    output_table.ttype{ch} = ttype;
    output_table.tstat(ch) = statvalue;
    output_table.p(ch) = p;
    output_table.h(ch) = h;
    output_table.p_FDR(ch) = p;
    output_table.h_FDR(ch) = 0;
    output_table.v_lo{ch} = v_lo;
    output_table.v_hi{ch} = v_hi;
end

% sorting outcomes in descending order or significance (1-p)
output_table = sortrows(output_table, 'p', 'ascend');

end
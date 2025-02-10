function [ output_table ] = func_pairwise_RSN_stat(data_hc, data_sz, brlab, meas)

% Function compare spectral indices obtained from HC and SZ groups
% (independent design) on the level of resting-state networks.
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 02/06/2025

%% preallocation and analysis constants

nr = length(brlab); % # of resting-state networks

output_table = table(cell(nr,1), cell(nr,1), zeros(nr,1),...
    zeros(nr,1), zeros(nr,1),...
    zeros(nr,1), zeros(nr,1),...
    cell(nr,1), zeros(nr,1),...
    zeros(nr,1), zeros(nr,1),...
    cell(nr,1), cell(nr,1),...
    'VariableNames', {'meas','rsn', 'rsnID', 'E_hc','E_sz','p','h','ttype','tstat','p_FDR','h_FDR','v_hc','v_sz'});

%% statistical tests

for br = 1:nr % iterate over RSNs

    % select variables
    v_hc = data_hc(:,br);
    v_sz = data_sz(:,br);

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
    output_table.meas{br} = meas;
    output_table.rsn{br} = brlab{br};
    output_table.rsnID(br) = br;
    output_table.E_hc(br) = E_hc;
    output_table.E_sz(br) = E_sz;
    output_table.ttype{br} = ttype;
    output_table.tstat(br) = statvalue;
    output_table.p(br) = p;
    output_table.h(br) = h;
    output_table.p_FDR(br) = p;
    output_table.h_FDR(br) = 0;
    output_table.v_hc{br} = v_hc;
    output_table.v_sz{br} = v_sz;
end

% sorting outcomes in descending order or significance (1-p)
output_table = sortrows(output_table, 'p', 'ascend');

end
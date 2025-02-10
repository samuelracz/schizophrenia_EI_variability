function [ output_table ] = func_pairwise_FDR(input_table, FDR_alpha)

% Function to adjust outcomes to multiple comparisons using the False
% Discovery Rate (FDR) method of Benjamini & Hochberg (1995)
%
% Note: This function expects the input in tabular format where adjusted p 
% and h values are preallocated with zeros. Outcome remains sorted in order
% of significance.
%
% Reference: Benjamini & Hochberg. "Controlling the false discovery rate:
% a practical and powerful approach to multiple testing." Journal of the
% Royal statistical society: series B (Methodological) (1995).
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 02/06/2025

% number of comaprisons
N = size(input_table,1);

% sorting outcomes in decreasing order of significance (1-p)
output_table = sortrows(input_table, 'p', 'ascend');

% adjusting p-values
for n = 1:N
    output_table.p_FDR(n) = output_table.p(n)*N/n;
end

% identifying significant differences
for n = 1:N % iterate over comparisons
    if output_table.p_FDR(n) < FDR_alpha
        output_table.h_FDR(n) = 1;
    else % first time the adjusted p-value > 0, terminate
        break
    end
end

end
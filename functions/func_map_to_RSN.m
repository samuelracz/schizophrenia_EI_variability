function [ M_out ] = func_map_to_RSN(M)

% Function to aggregate channels to reflect brain resting-state networks 
% (RSNs) as defined in Yeo et al. (2011). Channel allocation to the various
% RSN was based on Giacometti et al. (2014) and previously employed in Racz
% et al. (2019) and (2021).
%
% References:
% Yeo et al. "The organization of the human cerebral cortex estimated by 
% intrinsic functional connectivity." Journal of neurophysiology (2011).
%
% Giacometti et al. "Algorithm to find high density EEG scalp coordinates 
% and analysis of their correspondence to structural and functional regions
% of the brain." Journal of neuroscience methods (2014).
%
% Racz et al. "Multifractal and entropy analysis of resting-state
% electroencephalography reveals spatial organization in local dynamic
% functional connectivity." Scientific Reports (2019).
%
% Racz et al. "Separating scale‐free and oscillatory components of neural
% activity in schizophrenia." Brain and behavior (2021).
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 02/06/2025

% channel allocations:
VN = [37, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55]; % # of channels: 12
SM = [22, 23, 24, 25, 26, 27, 28, 32, 33, 34]; % of channels: 10
DA = [31, 35, 38, 39, 40, 41, 42, 43, 44]; % of channels: 9
VAL = [6, 12, 13, 14, 16, 17, 19, 20, 21, 29, 30, 36]; % of channels: 12
FP = [7, 8, 10, 11, 15, 18]; % of channels: 6
DMN = [1, 2, 3, 4, 5, 9]; % of channels: 6

labels = {'VN','SM','DA','VAL','FP','DMN'};
regs = {VN, SM, DA, VAL, FP, DMN};

N = length(regs);

M_out = struct('lab', {regs}, 'xyz', zeros(N,3));

for n = 1:N
    xyz_tmp = M.xyz(regs{n},:);
    M_out.xyz(n,:) = mean(xyz_tmp,1);
end

M_out.labels = labels;

end
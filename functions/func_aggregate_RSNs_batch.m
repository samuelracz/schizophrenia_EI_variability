function [ data_out ] = func_aggregate_RSNs_batch(data_in, Mr)

% Function to aggregate EEG channel indices over RSNs for multiple sliding
% windows
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 02/06/2025

regs = Mr.lab; %...........RSN names
T = size(data_in,1); %.....# of sliding windows
N = length(regs); %........# of RSNs

data_out = zeros(T,N);

for t = 1:N
    data_out(:,t) = mean(data_in(:,regs{t}),2);
end

end
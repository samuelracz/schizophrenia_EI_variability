function [ data_out, data_out_rsn ] = func_average_RSN_results_irasa(input_struct, chlab, Mr)

% Function to compute the average power spectra and spectral indices over
% sliding windows. Note that IRASA does not guarantee that the isolated
% fractal power is lower than the raw (mixed) spectral power and therefore
% the isolated oscillatory power (mixd-frac) can be negative, preventing
% log-transform. As we were primarily interested in log-transformed power,
% we insted computed the proportion of oscillatory power to total spectral
% power (i.e., log(osci)=log(mixd)-log(frac)).
%
% Note that variance over time in spectral power is also computed after
% log transformation.
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 02/06/2025


%% channel-wise averaging
data_out = struct(...
    'freq', [],...
    'srate', [],...
    'mixd', [],...
    'frac', [],...
    'osci', [],...
    'mixdvar', [],...
    'fracvar', [],...
    'oscivar', [],...
    'Beta_bb', [],...
    'Cons_bb', [],...
    'Plaw_bb', [],...
    'Freq_bb', [],...
    'Beta_lo', [],...
    'Cons_lo', [],...
    'Plaw_lo', [],...
    'Freq_lo', [],...
    'Beta_hi', [],...
    'Cons_hi', [],...
    'Plaw_hi', [],...
    'Freq_hi', []);

data_in = input_struct([input_struct.flag_af]==0);

%% analysis constants
N = length(data_in); % # of subjects
nch = length(chlab); % # of channels

%% preallocations
mat_mixd = [];
mat_frac = [];
mat_osci = [];
mat_beta_bb = zeros(nch,N);
mat_cons_bb = zeros(nch,N);
mat_beta_lo = zeros(nch,N);
mat_cons_lo = zeros(nch,N);
mat_beta_hi = zeros(nch,N);
mat_cons_hi = zeros(nch,N);
mat_plaw_bb = [];
mat_plaw_lo = [];
mat_plaw_hi = [];

%% average over sliding windows

% iterate over sliding windows
for n = 1:N
    tmp = data_in(n).plaw;
    
    % preallocate
    if isempty(mat_mixd) && isempty(mat_frac) && isempty(mat_osci)
        mat_mixd = zeros([size(tmp.mixd),N]);
        mat_frac = zeros([size(tmp.frac),N]);
        mat_osci = zeros([size(tmp.osci),N]);
    end

    % preallocate
    if isempty(mat_plaw_bb) && isempty(mat_plaw_lo) && isempty(mat_plaw_hi)
        mat_plaw_bb = zeros([size(tmp.Plaw_bb),N]);
        mat_plaw_lo = zeros([size(tmp.Plaw_lo),N]);
        mat_plaw_hi = zeros([size(tmp.Plaw_hi),N]);
    end

    % store data
    mat_mixd(:,:,n) = tmp.mixd;
    mat_frac(:,:,n) = tmp.frac;
    mat_osci(:,:,n) = exp(log(tmp.mixd) - log(tmp.frac));
    mat_beta_bb(:,n) = tmp.Beta_bb;
    mat_cons_bb(:,n) = tmp.Cons_bb;
    mat_plaw_bb(:,:,n) = tmp.Plaw_bb;
    mat_beta_lo(:,n) = tmp.Beta_lo;
    mat_cons_lo(:,n) = tmp.Cons_lo;
    mat_plaw_lo(:,:,n) = tmp.Plaw_lo;
    mat_beta_hi(:,n) = tmp.Beta_hi;
    mat_cons_hi(:,n) = tmp.Cons_hi;
    mat_plaw_hi(:,:,n) = tmp.Plaw_hi;
end

% average over sliding windows
data_out.freq = tmp.freq;
data_out.srate = tmp.srate;
data_out.mixd = mean(mat_mixd,3);
data_out.frac = mean(mat_frac,3);
data_out.osci = mean(mat_osci,3);
data_out.mixdvar = exp(std(log(mat_mixd),[],3));
data_out.fracvar = exp(std(log(mat_frac),[],3));
data_out.oscivar = exp(std(log(mat_osci),[],3));
data_out.Beta_bb = mean(mat_beta_bb,2);
data_out.Cons_bb = mean(mat_cons_bb,2);
data_out.Plaw_bb = mean(mat_plaw_bb,3);
data_out.Freq_bb = tmp.Freq_bb;
data_out.Beta_lo = mean(mat_beta_lo,2);
data_out.Cons_lo = mean(mat_cons_lo,2);
data_out.Plaw_lo = mean(mat_plaw_lo,3);
data_out.Freq_lo = tmp.Freq_lo;
data_out.Beta_hi = mean(mat_beta_hi,2);
data_out.Cons_hi = mean(mat_cons_hi,2);
data_out.Plaw_hi = mean(mat_plaw_hi,3);
data_out.Freq_hi = tmp.Freq_hi;

%% resting-state network-level averaging
data_out_rsn = struct(...
    'freq', [],...
    'srate', [],...
    'mixd', [],...
    'frac', [],...
    'osci', [],...
    'mixdvar', [],...
    'fracvar', [],...
    'oscivar', [],...
    'Beta_bb', [],...
    'Cons_bb', [],...
    'Plaw_bb', [],...
    'Freq_bb', [],...
    'Beta_lo', [],...
    'Cons_lo', [],...
    'Plaw_lo', [],...
    'Freq_lo', [],...
    'Beta_hi', [],...
    'Cons_hi', [],...
    'Plaw_hi', [],...
    'Freq_hi', []);

%% analysis constants
rsnlab = Mr.lab; %........RSN labels
nb = length(rsnlab); %....# of RSNs

%% preallocations
mat_mixd = [];
mat_frac = [];
mat_osci = [];

mat_mixdvar = [];
mat_fracvar = [];
mat_oscivar = [];

mat_beta_bb = zeros(nb,1);
mat_cons_bb = zeros(nb,1);
mat_beta_lo = zeros(nb,1);
mat_cons_lo = zeros(nb,1);
mat_beta_hi = zeros(nb,1);
mat_cons_hi = zeros(nb,1);

mat_plaw_bb = [];
mat_plaw_lo = [];
mat_plaw_hi = [];

%% average over RSNs

% iterate over RSNs
for br = 1:length(rsnlab)
    % preallocate
    if isempty(mat_mixd) && isempty(mat_frac) && isempty(mat_osci)
        mat_mixd = zeros([size(data_out.mixd,1),nb]);
        mat_frac = zeros([size(data_out.frac,1),nb]);
        mat_osci = zeros([size(data_out.osci,1),nb]);
        mat_mixdvar = zeros([size(data_out.mixdvar,1),nb]);
        mat_fracvar = zeros([size(data_out.fracvar,1),nb]);
        mat_oscivar = zeros([size(data_out.oscivar,1),nb]);
    end

    % preallocate
    if isempty(mat_plaw_bb) && isempty(mat_plaw_lo) && isempty(mat_plaw_hi)
        mat_plaw_bb = zeros([size(data_out.Plaw_bb,1),nb]);
        mat_plaw_lo = zeros([size(data_out.Plaw_lo,1),nb]);
        mat_plaw_hi = zeros([size(data_out.Plaw_hi,1),nb]);
    end

    % aggregate channel indices over RSNs
    mat_mixd(:,br) = (mean((data_out.mixd(:,rsnlab{br})),2));
    mat_frac(:,br) = (mean((data_out.frac(:,rsnlab{br})),2));
    mat_osci(:,br) = (mean((data_out.mixd(:,rsnlab{br}))-log(data_out.frac(:,rsnlab{br})),2));

    mat_mixdvar(:,br) = (mean((data_out.mixdvar(:,rsnlab{br})),2));
    mat_fracvar(:,br) = (mean((data_out.fracvar(:,rsnlab{br})),2));
    mat_oscivar(:,br) = (mean((data_out.oscivar(:,rsnlab{br})),2));

    mat_beta_bb(br) = mean(data_out.Beta_bb(rsnlab{br}));
    mat_cons_bb(br) = mean(data_out.Cons_bb(rsnlab{br}));
    mat_plaw_bb(:,br) = mean(data_out.Plaw_bb(:,rsnlab{br}),2);
    mat_beta_lo(br) = mean(data_out.Beta_lo(rsnlab{br}));
    mat_cons_lo(br) = mean(data_out.Cons_lo(rsnlab{br}));
    mat_plaw_lo(:,br) = mean(data_out.Plaw_lo(:,rsnlab{br}),2);
    mat_beta_hi(br) = mean(data_out.Beta_hi(rsnlab{br}));
    mat_cons_hi(br) = mean(data_out.Cons_hi(rsnlab{br}));
    mat_plaw_hi(:,br) = mean(data_out.Plaw_hi(:,rsnlab{br}),2);
    
end

% store data
data_out_rsn.freq = data_out.freq;
data_out_rsn.srate = data_out.srate;
data_out_rsn.mixd = mat_mixd;
data_out_rsn.frac = mat_frac;
data_out_rsn.osci = mat_osci;
data_out_rsn.mixdvar = mat_mixdvar;
data_out_rsn.fracvar = mat_fracvar;
data_out_rsn.oscivar = mat_oscivar;
data_out_rsn.Beta_bb = mat_beta_bb;
data_out_rsn.Cons_bb = mat_cons_bb;
data_out_rsn.Plaw_bb = mat_plaw_bb;
data_out_rsn.Freq_bb = data_out.Freq_bb;
data_out_rsn.Beta_lo = mat_beta_lo;
data_out_rsn.Cons_lo = mat_cons_lo;
data_out_rsn.Plaw_lo = mat_plaw_lo;
data_out_rsn.Freq_lo = data_out.Freq_lo;
data_out_rsn.Beta_hi = mat_beta_hi;
data_out_rsn.Cons_hi = mat_cons_hi;
data_out_rsn.Plaw_hi = mat_plaw_hi;
data_out_rsn.Freq_hi = data_out.Freq_hi;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis script to compute spectral slopes in a dynamic manner using
% Irregular Resampling Auto-Spectral Analysis (IRASA) by Wen & Liu (2016).
% Analysis is performed on the same datasets evaluated in Racz et al.
% (2025); for details on pre-processing, see script_00_preproc_ec_mara.m
%
% NOTE: this script loads the pre-processed data from MATLAB workspaces,
% however corresponding .edf files are also provided.
%
% Reference for the IRASA method:
% Wen H. and Liu Z. (2015) Separating Fractal and Oscillatory Components in
% the Power Spectrum of Neurophysiological Signals
% 
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 02/06/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('functions'))
addpath('miscellaneous')

OD = cd;
fnames_ec = dir('data_ec_mara_avg_256Hz_mat/*.mat');

%% analysis constants
nch = 55; %...................number of channels
fs = 256; %...................frequency of sampling
ns = length(fnames_ec); %.....total number of participants
frange_bb = [1 45]; %.........broadband frequency range
frange_lo = [1, 4]; %.........low-frequency range
frange_hi = [20 45];%.........high-frequency range
hset = 1.1:0.1:2.6; %.........set of resampling factors

W = 8; %......................sliding window analysis window size (s)
wstep = 2; %..................sliding window analysis step size (s)
nseg = (30-W)/wstep+1; %......total number of windows

%% analysis
for subj = 1:ns
    tic
    % load data
    load(['data_ec_mara_avg_256Hz_mat/' fnames_ec(subj).name])
    eeg_ec = data_ec(1:nch,:); % remove trigger channel
    
    chlab = cell(1,nch);
    for ch = 1:nch
        chlab{ch} = hdr.label{ch};
    end

    % IRASA from full signal (30s) - only for sanity check
    spec_ec = amri_sig_fractal(eeg_ec',fs,'hset',hset);
    plaw_ec = func_IRASA_plawfit_multi(spec_ec, frange_bb, frange_lo, frange_hi);
    
    % reshape data
    eeg_seg_ec = zeros(nch,W*fs,nseg);
    for seg = 1:nseg
        eeg_seg_ec(:,:,seg) = eeg_ec(:,(seg-1)*wstep*fs+1:(seg-1)*wstep*fs+W*fs);
    end
    
    % IRASA from segmented data (8s segments, 75% overlap)
    results_seg_ec = struct('ID',[],'spec',[],'plaw',[],'flag_af',[]);
    for seg = 1:nseg
        spec_seg = amri_sig_fractal(squeeze(eeg_seg_ec(:,:,seg))',fs,'hset',hset);
        plaw_seg = func_IRASA_plawfit_multi(spec_seg,frange_bb,frange_lo,frange_hi);

        % store output
        results_seg_ec(seg).ID = seg;
        results_seg_ec(seg).spec = spec_seg;
        results_seg_ec(seg).plaw = plaw_seg;
        results_seg_ec(seg).flag_af = 0;
    end

    % saving results
    if ~exist('results_IRASA_ec', 'dir')
       mkdir('results_IRASA_ec')
    end
    subj_saved = subj;
    hdr_saved = hdr;
    fname_ec = [fnames_ec(subj).name(1:end-8) '_IRASA_res.mat'];
    save(['results_IRASA_ec/' fname_ec], 'spec_ec', 'plaw_ec', 'results_seg_ec', 'hdr_saved', 'subj_saved');
    disp([fname_ec, ' || time elapsed: ', num2str(toc)])
end
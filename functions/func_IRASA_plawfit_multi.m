function [ plaw ] = func_IRASA_plawfit_multi(spec, frange_bb, frange_lo, frange_hi)

% This function is derived from (and is utilizing) amri_sig_plawfit.m from 
% the IRASA toolbox by Wen & Liu (2016), tailored to fit the spectral slope
% in multiple frequency ranges (broadband, low- and high-frequency).
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 02/06/2025
%
% Reference for IRASA:
% Wen H. and Liu Z. (2015) Separating Fractal and Oscillatory Components in
% the Power Spectrum of Neurophysiological Signals

plaw = spec;

nch = size(spec.mixd,2);

beta_bb = zeros(nch,1);
cons_bb = zeros(nch,1);
beta_lo = zeros(nch,1);
cons_lo = zeros(nch,1);
beta_hi = zeros(nch,1);
cons_hi = zeros(nch,1);
pl_bb = [];
pl_lo = [];
pl_hi = [];

for ch = 1:nch
    spec_tmp = struct(...
        'freq', spec.freq,...
        'srate', spec.srate,...
        'mixd', spec.mixd(:,ch),...
        'frac', spec.frac(:,ch),...
        'osci', spec.osci(:,ch));

    % broadband
    plaw_bb = amri_sig_plawfit(spec_tmp, frange_bb);
    beta_bb(ch) = plaw_bb.Beta;
    cons_bb(ch) = plaw_bb.Cons;

    if isempty(pl_bb)
        pl_bb = zeros(length(plaw_bb.Plaw), nch);
    end
    pl_bb(:,ch) = plaw_bb.Plaw;
    fr_bb = plaw_bb.Freq;

    % low-frequency
    plaw_lo = amri_sig_plawfit(spec_tmp, frange_lo);
    beta_lo(ch) = plaw_lo.Beta;
    cons_lo(ch) = plaw_lo.Cons;

    if isempty(pl_lo)
        pl_lo = zeros(length(plaw_lo.Plaw),nch);
    end
    pl_lo(:,ch) = plaw_lo.Plaw;
    fr_lo = plaw_lo.Freq;

    % high-frequency
    plaw_hi = amri_sig_plawfit(spec_tmp, frange_hi);
    beta_hi(ch) = plaw_hi.Beta;
    cons_hi(ch) = plaw_hi.Cons;

    if isempty(pl_hi)
        pl_hi = zeros(length(plaw_hi.Plaw),nch);
    end
    pl_hi(:,ch) = plaw_hi.Plaw;
    fr_hi = plaw_hi.Freq;
end

% broadband
plaw.Beta_bb = beta_bb;
plaw.Cons_bb = cons_bb;
plaw.Plaw_bb = pl_bb;
plaw.Freq_bb = fr_bb;

% low-frequency
plaw.Beta_lo=  beta_lo;
plaw.Cons_lo = cons_lo;
plaw.Plaw_lo = pl_lo;
plaw.Freq_lo = fr_lo;

% high-frequency
plaw.Beta_hi = beta_hi;
plaw.Cons_hi = cons_hi;
plaw.Plaw_hi = pl_hi;
plaw.Freq_hi = fr_hi;

% Compares the following filters for combined noise reduction (NR) and
% acoustic echo cancellation (AEC):
% *) MWFext: An extended multichanncel Wiener filter [1,2]
% *) GEIC: A generalised acoustic echo and interference canceller [3]
%
% v1.0
% LICENSE: This software is distributed under the terms of the MIT license (See LICENSE.md).
% AUTHOR:  Arnout Roebben
% CONTACT: arnout.roebben@esat.kuleuven.be
% CITE: A. Roebben, T. van Waterschoot, and M. Moonen, "A comparative 
% analysis of generalised echo and interference cancelling and extended 
% multichannel Wiener filtering for combined noise reduction and acoustic
% echo cancellation, Accepted for publication in 2025 IEEE
% International Conference on Acoustics, Speech and Signal Processing
% (ICASSP), Hyderabad, India, Apr. 2025.
% and
% A. Roebben, “Github repository: A Comparative analysis of
% generalised echo and interference cancelling and extended
% multichannel Wiener filtering for combined noise reduction
% and acoustic echo cancellation,” https://https://github.com/Arnout-
% Roebben/AECNR_GEIC_vs_MWFext, 2025.
%
% A preprint is available at
% A. Roebben, T. van Waterschoot, and M. Moonen, "A comparative 
% analysis of generalised echo and interference cancelling and extended 
% multichannel Wiener filtering for combined noise reduction and acoustic
% echo cancellation, 2025, arxiv:2503.03593.
%
% [1] Geert Rombouts and Marc Moonen, “An integrated approach to acoustic 
%     noise and echo cancellation,” Signal Processing, vol. 85, no. 4, 
%     pp. 849–871, Apr. 2005.
% [2] Santiago Ruiz, Toon van Waterschoot, and Marc Moonen, “Distributed 
%     Combined Acoustic Echo Cancellation and Noise Reduction in Wireless 
%     Acoustic Sensor and Actuator Networks,” IEEE/ACM Transactions on 
%     Audio, Speech, and Language Processing, vol. 30, pp. 534–547, 2022.
% [3] W. Herbordt, W. Kellermann, and S. Nakamura, “Joint optimization of 
%     LCMV beamforming and acoustic echo cancellation,” in 2004 12th 
%     European Signal Processing Conference (EUSIPCO), Vienna, Austria, 
%     Sept. 2004, pp. 2003–2006.
% [4] S. Markovich-Golan and S. Gannot, “Performance analysis of the
%     covariance subtraction method for relative transfer function estimation
%     and comparison to the covariance whitening method,” in 2015 IEEE
%     International Conference on Acoustics, Speech and Signal Processing
%     (ICASSP), Apr. 2015, pp. 544–548.

%% Check that current folder corresponds to AECNR_GEIC_vs_MWFext
[~,curr_dir] = fileparts(pwd);
if ~strcmp(curr_dir,"AECNR_GEIC_vs_MWFext")
    error('Current fulder must correspond to ''AECNR_GEIC_vs_MWFext''!')
end

%% Cleanup
clear; close all; clc;
rng(2,"twister"); % Fix random number generator
addpath(genpath('.')); % Add subfolders of directory to path 

%% Load audio
% Audio data consists of the following components
% -m    TXM     M-microphone microphone signal of length T samples.
%               m=s+n+e.
% -s    TXM     M-microphone desired speech signal of length T samples.
% -n    TXM     M-microphone near-end room noise signal of length T samples.
% -e    TXM     M-microphone far-end room speech component in the echo 
%               signal of length T samples.
% -l    TXL     L-loudspeaker loudspeaker signal of length T samples. 
% For the desired speech and far-end room speech component in the
% loudspeakers, sentences from the Voice Cloning Toolkit (VCTK) corpus were
% used [1], as also available at [2].
% [1] C. Veaux, J. Yamagishi, K. MacDonald, "CSTR VCTK Corpus:
% English Multi-speaker Corpus for CSTR Voice Cloning Toolkit," http:
% //homepages.inf.ed.ac.uk/jyamagis/page3/page58/page58.html, 2016.
% [2] Dietzen, T., Ali, R., Taseska, M., van Waterschoot, T.: "Data
% Repository for MYRiAD: A Multi-Array Room Acoustic Database,"
% https://zenodo.org/record/7389996, 2023.
load('.\Audio\sig.mat');
% Load desired speech source-to-mic impulse response
load('.\Audio\impulse.mat'); 

%% Processing parameters
p = struct();

% General 
p.fs = fs; % Sampling rate
p.ref = 1; % Reference microphone
p.M = size(sig.m,2); % Amount of microphones
p.L = size(sig.l,2); % Amount of loudspeakers

% Desired speech source-to-mic impulse response
p.hs = hs;

% Frequency transform (See also WOLA_analysis.m and WOLA_synthesis.m)
p.N = 512; % Discrete Fourier transform (DFT) size N
p.shift = p.N/2; % Frame shift for weighted overlap add (WOLA)
p.win = sqrt(hann(p.N,'periodic')); % Window 

% Voice activity detection (VAD) (See also VAD.m)
p.sensitivity = 1e-5; % Sensitivity of VAD

% Parameters (See also process_GEIC.m and process_MWFext.m)
p.rank_s = 1; % Requested rank of desired speech correlation matrix 
p.lambda = 0.995; % Smoothing factor when updating the matrices

% Metrics (See align_proc_unproc.m)
% Start time [s], after which the data is used to compute the metrics.
p.T_start = 5; 

%% Process
% MWFext [1,2] using a generalised eigenvalue decomposition (GEVD) approximation
% for the desired speech correlation matrix [2]
MWFext = process(sig,p,'MWFext');
% GEIC [3] using the ground truth relative transfer functions 
GEIC = process(sig,p,'GEIC');
% GEIC [3] using the estimated relative transfer functions using the
% covariance whitening (CW) method
GEIC_GEVD = process(sig,p,'GEIC_GEVD');

%% Metrics
[metrics_ref,MWFext.metrics] = calculateMetrics(sig,MWFext,p);  
[~,GEIC.metrics] = calculateMetrics(sig,GEIC,p);
[~,GEIC_GEVD.metrics] = calculateMetrics(sig,GEIC_GEVD,p);

%% Visualisation
% Metrics
fprintf('MWFext:\n')
fprintf('\t SNR improvement: %f\n',MWFext.metrics.snr - metrics_ref.snr);
fprintf('\t SD: %f\n',MWFext.metrics.sd);
fprintf('\t ERLE: %f\n\n',MWFext.metrics.erle);

fprintf('GEIC:\n')
fprintf('\t SNR improvement: %f\n',GEIC.metrics.snr - metrics_ref.snr);
fprintf('\t SD: %f\n',GEIC.metrics.sd);
fprintf('\t ERLE: %f\n\n',GEIC.metrics.erle);

fprintf('GEIC-GEVD:\n')
fprintf('\t SNR improvement: %f\n',GEIC_GEVD.metrics.snr - metrics_ref.snr);
fprintf('\t SD: %f\n',GEIC_GEVD.metrics.sd);
fprintf('\t ERLE: %f\n',GEIC_GEVD.metrics.erle);

% Signals
figure; hold on
t = tiledlayout(3,1);
ax = nexttile; hold on; plot(sig.s(:,p.ref)); plot(MWFext.time.s(:,p.ref)); 
title('Desired speech'); xlabel('Time [Samples]'); ylabel('Amplitude [Arb. unit]');
nexttile; hold on; plot(sig.n(:,p.ref)); plot(MWFext.time.n(:,p.ref)); 
title('Noise'); xlabel('Time [Samples]'); ylabel('Amplitude [Arb. unit]');
nexttile; hold on; plot(sig.e(:,p.ref)); plot(MWFext.time.e(:,p.ref)); 
title('Echo'); xlabel('Time [Samples]'); ylabel('Amplitude [Arb. unit]');
title(t,'MWFext');
lg  = legend(ax,["Input" "Output"],'Orientation','Horizontal'); 
lg.Layout.Tile = 'South'; hold off;

figure; hold on
t = tiledlayout(3,1);
ax = nexttile; hold on; plot(sig.s(:,p.ref)); plot(GEIC.time.s(:,p.ref)); 
title('Desired speech'); xlabel('Time [Samples]'); ylabel('Amplitude [Arb. unit]');
nexttile; hold on; plot(sig.n(:,p.ref)); plot(GEIC.time.n(:,p.ref)); 
title('Noise'); xlabel('Time [Samples]'); ylabel('Amplitude [Arb. unit]');
nexttile; hold on; plot(sig.e(:,p.ref)); plot(GEIC.time.e(:,p.ref)); 
title('Echo'); xlabel('Time [Samples]'); ylabel('Amplitude [Arb. unit]');
title(t,'GEIC');
lg  = legend(ax,["Input" "Output"],'Orientation','Horizontal'); 
lg.Layout.Tile = 'South'; hold off;

figure; hold on
t = tiledlayout(3,1);
ax = nexttile; hold on; plot(sig.s(:,p.ref)); plot(GEIC_GEVD.time.s(:,p.ref)); 
title('Desired speech'); xlabel('Time [Samples]'); ylabel('Amplitude [Arb. unit]');
nexttile; hold on; plot(sig.n(:,p.ref)); plot(GEIC_GEVD.time.n(:,p.ref)); 
title('Noise'); xlabel('Time [Samples]'); ylabel('Amplitude [Arb. unit]');
nexttile; hold on; plot(sig.e(:,p.ref)); plot(GEIC_GEVD.time.e(:,p.ref)); 
title('Echo'); xlabel('Time [Samples]'); ylabel('Amplitude [Arb. unit]');
title(t,'GEIC-GEVD');
lg  = legend(ax,["Input" "Output"],'Orientation','Horizontal'); 
lg.Layout.Tile = 'South'; hold off;
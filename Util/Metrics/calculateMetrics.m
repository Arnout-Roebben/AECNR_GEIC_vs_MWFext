function [metrics_ref,metrics_processed] = calculateMetrics(sig,processed,param)
% Computes the metrics to evaluate the signals in sig.
%
% INPUT: 
% sig           Struct      Struct containing the following input signals:
% -m            TXM         M-microphone microphone signal of length T samples.
%                           m=s+n+e.
% -s            TXM         M-microphone desired speech signal of length T samples.
% -n            TXM         M-microphone noise signal of length T samples.
% -e            TXM         M-microphone echo signal of length T samples.
%
% OUTPUT:
% res           Struct      Struct containing the computed metrics:
% -snr          1X1         Signal to noise ratio.
% -sd           1X1         Speech distortion
% -erle         1X1         Echo return loss enhancement
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

%% Initialisation
% Throw error if start time exceeds signal length
if size(sig.s,1)/param.fs < param.T_start
    error('The signals end before the supplied ''T_start''!');
end

% Adjust starting and end time of...
% ... Input signals
sm = sig.s(floor(param.T_start*param.fs):end,:);
mm = sig.m(floor(param.T_start*param.fs):end,:);
nm = sig.n(floor(param.T_start*param.fs):end,:);
em = sig.e(floor(param.T_start*param.fs):end,:);

% ... Output signals
s = processed.time.s(floor(param.T_start*param.fs):end,:);
m = processed.time.m(floor(param.T_start*param.fs):end,:);
n = processed.time.n(floor(param.T_start*param.fs):end,:);
e = processed.time.e(floor(param.T_start*param.fs):end,:);

% Select frames where speech is active of...
VADs_time = abs(sm(:,param.ref)) > std(sm(:,param.ref))*param.sensitivity;
% ... Input signals
sm = sm(VADs_time,param.ref);
mm = mm(VADs_time,param.ref);
nm = nm(VADs_time,param.ref);
em = em(VADs_time,param.ref);

% ... Output signals
s = s(VADs_time,:);
m = m(VADs_time,:);
n = n(VADs_time,:);
e = e(VADs_time,:);

%% Reference metric: without processing
metrics_ref = struct();
metrics_ref.snr = SNR(sm(:,param.ref),nm(:,param.ref));

%% Metrics: after processing
metrics_processed = struct();
metrics_processed.snr = SNR(s,n);
metrics_processed.sd = SD(sm,s);
metrics_processed.erle = ERLE(em, e);
end
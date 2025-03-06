function processed = process(sig,p,mode)
% Processing using an extended multichannel Wiener filter (MWFext) [1,2], 
% a generalised echo and interference canceller [3].
%
% INPUT: 
% sig           Struct      Struct containing the following input signals:
% -m            TXM         M-microphone microphone signal of length T samples.
%                           m=s+n+es+en.
% -s            TXM         M-microphone desired speech signal of length T samples.
% -n            TXM         M-microphone near-end room noise signal of length T samples.
% -e            TXM         M-microphone far-end room echo signal of length T samples.
% -l            TXL         L-loudspeaker loudspeaker signal of length T samples. 
% p             Struct      Struct containing the following parameters:
% -ref          1X1         Reference microphone.
% -rank_s       1X1         [Optional] Rank to be used in the 'GEVD' procedure 
%                           for the desired speech correlation matrix. See
%                           updateDifferenceCorrelation.m.
% -sensitivity  String      Sensitivity of the standard deviation in the 
%                           voice acitivity detector (VAD) formula, 
%                           see VAD.m.
% -fs           1X1         Sampling rate [Hz].
% -M            1X1         Number of microphones.
% -L            1X1         Number of loudspeakers.
% -N            1X1         Discrete Fourier transform (DFT) size. 
%                           See WOLA_analysis.m  and WOLA_synthesis.m
% -win          NX1         Window. See WOLA_analysis.m and WOLA_synthesis.m
% -shift        1X1         Frame shift. See WOLA_analysis.m and WOLA_synthesis.m
% -lambda       1X1         Exponential smoothing factor for correlaton
%                           matrix averaging, signifying the weight associated 
%                           to the previous estimate.
% -hs           FXM         M desired speech source-to-microphone transfer
%                           functions of length F.
%
% OUTPUT:         
% processed     Struct      Struct containing the processed signals, which is either:
% -MWFext       Struct      Struct containing the processed signals after
%                           MWFext.
% -GEIC         Struct      Struct containing the processed signals after
%                           GEIC using ground truth relative transfer functions.
% -GEIC_GEVD    Struct      Struct containing the processed signals after
%                           GEIC using the covariance whitening method [4] to
%                           estimate the relative transfer functions.
% --time        Struct      Struct containing the output in time domain for
%                           the reference microphone.
%                          
% ---m          TX1         See INPUT.
% ---s          TX1         See INPUT.
% ---n          TX1         See INPUT.
% ---e          TX1         See INPUT.
% 
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

%% Check mode argument
if ~strcmp(mode,'MWFext') && ~strcmp(mode,'GEIC')  && ~strcmp(mode,'GEIC_GEVD')
        error('Supplied mode not supported!');
end

%% Parse signals input
m = sig.m; 
s = sig.s;
n = sig.n;
l = sig.l;
e = sig.e;

%% STFT transformation: M X number of frames K X (N/2+1)
m_f = WOLA_analysis(m,p.win,p.N,p.shift);
s_f = WOLA_analysis(s,p.win,p.N,p.shift);
n_f = WOLA_analysis(n,p.win,p.N,p.shift);
l_f = WOLA_analysis(l,p.win,p.N,p.shift);
e_f = WOLA_analysis(e,p.win,p.N,p.shift);

%% Voice activity detection (VAD) calculation
VADs = VAD(s_f,p.sensitivity,p.ref); % VAD speech

%% Prepare signals
signals = struct();
signals.m_f = m_f; signals.s_f = s_f; signals.n_f = n_f; 
signals.e_f = e_f; signals.l_f = l_f; 

%% MWFext 
if strcmp(mode,'MWFext')
    % Parameters for processing
    param = struct();
    param.rank = p.rank_s; param.VADs = VADs; 
    param.lambda = p.lambda;

    % Processing
    processed = process_MWFext(signals,param);

    % Conversion to time domain
    processed.time.m = WOLA_synthesis(processed.m_f(p.ref,:,:),p.win,p.N,p.shift);
    processed.time.s = WOLA_synthesis(processed.s_f(p.ref,:,:),p.win,p.N,p.shift);
    processed.time.n = WOLA_synthesis(processed.n_f(p.ref,:,:),p.win,p.N,p.shift);
    processed.time.e = WOLA_synthesis(processed.e_f(p.ref,:,:),p.win,p.N,p.shift);

    % Remove STFT signals from output struct
    processed = rmfield(processed,'m_f'); processed = rmfield(processed,'s_f'); 
    processed = rmfield(processed,'n_f'); processed = rmfield(processed,'e_f'); 
end

%% GEIC
if strcmp(mode,'GEIC') 
    % Parameters for processing
    param = struct();
    param.VADs = VADs;
    param.lambda = p.lambda;

    % Compute Relative-transfer-functions
    param.Hs = fft(p.hs,p.N);
    param.Hs = permute(param.Hs(1:p.N/2+1,:,:),[2 3 1]);
    for n=1:p.N/2+1 
        param.Hs(:,:,n) = param.Hs(:,:,n)./param.Hs(p.ref,:,n); 
    end

    % Processing
    processed = process_GEIC(signals,param);

    % Conversion to time domain
    processed.time.m = WOLA_synthesis(processed.m_f,p.win,p.N,p.shift);
    processed.time.s = WOLA_synthesis(processed.s_f,p.win,p.N,p.shift);
    processed.time.n = WOLA_synthesis(processed.n_f,p.win,p.N,p.shift);
    processed.time.e = WOLA_synthesis(processed.e_f,p.win,p.N,p.shift); 

    % Remove STFT signals from output struct
    processed = rmfield(processed,'m_f'); processed = rmfield(processed,'s_f'); 
    processed = rmfield(processed,'n_f'); processed = rmfield(processed,'e_f');     
end

%% GEIC-GEVD
if strcmp(mode,'GEIC_GEVD') 
    % Parameters for processing
    param = struct();
    param.VADs = VADs;
    param.lambda = p.lambda;
    param.ref = p.ref;
    
    % Processing
    processed = process_GEIC_GEVD(signals,param);

    % Conversion to time domain
    processed.time.m = WOLA_synthesis(processed.m_f,p.win,p.N,p.shift);
    processed.time.s = WOLA_synthesis(processed.s_f,p.win,p.N,p.shift);
    processed.time.n = WOLA_synthesis(processed.n_f,p.win,p.N,p.shift);
    processed.time.e = WOLA_synthesis(processed.e_f,p.win,p.N,p.shift); 

    % Remove STFT signals from output struct
    processed = rmfield(processed,'m_f'); processed = rmfield(processed,'s_f'); 
    processed = rmfield(processed,'n_f'); processed = rmfield(processed,'e_f');     
end


end
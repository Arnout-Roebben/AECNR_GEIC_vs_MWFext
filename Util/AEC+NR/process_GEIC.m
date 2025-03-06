function GEIC = process_GEIC(signals,param)
% Processing using a generalised echo and interference canceller [1] using the ground truth
% relative transfer functions.
%
% INPUT:
% signals      Struct     STFT-transformed input signals:
% -m_f         MXKXN      M-microphone microphone signal with K frames and N frequency bins.
% -s_f         MXKXN      M-microphone desired speech signal with K frames and N frequency bins.
% -n_f         MXKXN      M-microphone noise signal with K frames and N frequency bins.
% -e_f         MXKXN      M-microphone echo signal with K frames and N frequency bins.
% -s_f         LXKXN      M-microphone loudspeaker signal with K frames and N frequency bins.
% -Hs          MX1XN      M STFT-transformed desired speech
%                         source-to-microphone transfer functions for N frequency bins.
% -lambda      1X1        Exponential smoothing factor for correlaton
%                         matrix averaging, signifying the weight associated 
%                         to the previous estimate.
% -VADs        KXN        Voice activity detector (VAD). A 1 denotes
%                         desired speech activity and a 0 desired speeech inacitvity.
%
% OUTPUT
% GEIC        Struct     STFT-transformed output signals:
% -m_f        MXKXN      See INPUT
% -s_f        MXKXN      See INPUT
% -n_f        MXKXN      See INPUT
% -e_f        MXKXN      See INPUT
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
% [1] W. Herbordt, W. Kellermann, and S. Nakamura, “Joint optimization of 
%     LCMV beamforming and acoustic echo cancellation,” in 2004 12th 
%     European Signal Processing Conference (EUSIPCO), Vienna, Austria, 
%     Sept. 2004, pp. 2003–2006.

%% Initialisation
M = size(signals.m_f,1); % Number of microphones
K = size(signals.m_f,2); % Number of frames
N = size(signals.m_f,3); % Number of channels
L = size(signals.l_f,1); % Number of loudspeakers

%% Construct blocking matrix
B = cell(N,1);
for n=1:N
    B{n} = null(param.Hs(:,:,n)')';
end

%% Construct quiescent beamformer
Fq = cell(N,1);
for n=1:N
    Fq{n} = param.Hs(:,:,n)/(param.Hs(:,:,n)'*param.Hs(:,:,n));
end

%% Adaptive filter
% Initialise correlation matrix
Rmm_f(:,:,:,1) = repmat(1e-6*eye(M+L),[1 1 N]);

for k = 1:K % Loop over frames
    for n=1:N % Loop over bins
        % Update correlation matrices
        if param.VADs(k,n)
            Rmm_f(:,:,n) = param.lambda*Rmm_f(:,:,n) + ...
                (1-param.lambda)*((cat(1,squeeze(signals.m_f(:,k,n)),...
                permute(signals.l_f(:,k,n),[1 3 2])))*(cat(1,squeeze(signals.m_f(:,k,n)),...
                permute(signals.l_f(:,k,n),[1 3 2])))');
        end
    end

    % Construct adaptive filter
    w = nan(M+L-1,N);
    for n=1:N
        w(:,n) =  pinv([B{n}*Rmm_f(1:M,1:M,n)*B{n}' B{n}*Rmm_f(1:M,M+1:end,n); Rmm_f(M+1:end,1:M,n)*B{n}' Rmm_f(M+1:end,M+1:end,n)])*[B{n}*Rmm_f(1:M,1:M,n)*Fq{n};Rmm_f(M+1:end,1:M,n)*Fq{n}];
    end

    % Apply filter
    % Fixed
    for n=1:N
        GEIC.m_f(:,k,n) = Fq{n}'*signals.m_f(:,k,n) - w(:,n)'*[B{n}*signals.m_f(:,k,n);signals.l_f(:,k,n)];
        GEIC.s_f(:,k,n) = Fq{n}'*signals.s_f(:,k,n) - w(:,n)'*[B{n}*signals.s_f(:,k,n);zeros(size(signals.l_f(:,k,n)))];
        GEIC.n_f(:,k,n) = Fq{n}'*signals.n_f(:,k,n) - w(:,n)'*[B{n}*signals.n_f(:,k,n);zeros(size(signals.l_f(:,k,n)))];
        GEIC.e_f(:,k,n) = Fq{n}'*signals.e_f(:,k,n) - w(:,n)'*[B{n}*signals.e_f(:,k,n);signals.l_f(:,k,n)];
    end
end

end
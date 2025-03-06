function GEIC_GEVD = process_GEIC_GEVD(signals,param)
% Processing using a generalised echo and interference canceller [1] using the 
% covariance whitening method [2] to estimate the relative transfer functions.
%
% INPUT:
% signals      Struct     STFT-transformed input signals:
% -m_f         MXKXN      M-microphone microphone signal with K frames and N frequency bins.
% -s_f         MXKXN      M-microphone desired speech signal with K frames and N frequency bins.
% -n_f         MXKXN      M-microphone noise signal with K frames and N frequency bins.
% -e_f         MXKXN      M-microphone echo signal with K frames and N frequency bins.
% -s_f         LXKXN      M-microphone loudspeaker signal with K frames and N frequency bins.
% -lambda      1X1        Exponential smoothing factor for correlaton
%                         matrix averaging, signifying the weight associated
%                         to the previous estimate.
% -VADs        KXN        Voice activity detector (VAD). A 1 denotes
%                         desired speech activity and a 0 desired speeech inacitvity.
% -ref         1X1        Reference microphone number
%
% OUTPUT
% GEIC_GEVD   Struct     STFT-transformed output signals:
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
% [2] S. Markovich-Golan and S. Gannot, “Performance analysis of the
%     covariance subtraction method for relative transfer function estimation
%     and comparison to the covariance whitening method,” in 2015 IEEE
%     International Conference on Acoustics, Speech and Signal Processing
%     (ICASSP), Apr. 2015, pp. 544–548.

%% Initialisation
M = size(signals.m_f,1); % Number of microphones
K = size(signals.m_f,2); % Number of frames
N = size(signals.m_f,3); % Number of channels
L = size(signals.l_f,1); % Number of loudspeakers

%% Collect correlation matrices
% Initialise correlation matrix when desired speech is active
Rmm1_f(:,:,:,1) = repmat(1e-6*eye(M+L),[1 1 N]);
% Initialise correlation matrix when desired speech is inactive
Rmm0_f(:,:,:,1) = repmat(1e-6*eye(M+L),[1 1 N]);
Rmm0_f(M+1:end,:,:,1) = 0;
Rmm0_f(:,M+1:end,:,1) = 0;

for k = 1:K % Loop over frames
    for n=1:N % Loop over bins
        % Update extended microphone signal whenever VADs(k,n)=1
        if param.VADs(k,n)
            Rmm1_f(:,:,n) = param.lambda*Rmm1_f(:,:,n) + ...
                (1-param.lambda)*((cat(1,squeeze(signals.m_f(:,k,n)),...
                permute(signals.l_f(:,k,n),[1 3 2])))*(cat(1,squeeze(signals.m_f(:,k,n)),...
                permute(signals.l_f(:,k,n),[1 3 2])))');
        else
            % Update extended microphone signal whenever VADs(k,n)=0
            Rmm0_f(:,:,n) = param.lambda*Rmm0_f(:,:,n) + ...
                (1-param.lambda)*((cat(1,squeeze(signals.m_f(:,k,n)),...
                permute(signals.l_f(:,k,n),[1 3 2])))*(cat(1,squeeze(signals.m_f(:,k,n)),...
                permute(signals.l_f(:,k,n),[1 3 2])))');
            Rmm0_f(M+1:end,:,:,n) = 0;
            Rmm0_f(:,M+1:end,:,n) = 0;
        end
    end

    % Construct blocking matrix and fixed beamformer
    Fq = cell(N,1);
    B = cell(N,1);
    for n=1:N
        % Get generalised eigenvectors
        [~,~,Q] = updateDifferenceCorrelation(Rmm1_f(:,:,n),Rmm0_f(:,:,n),1);
        % Get generalised eigenvectors with associated zeros to final L
        % elements
        [~,I] = sort(sum(abs(Q(M+1:end,:)),1),'ascend');
        I = sort(I(1:M),'ascend');
        % From eigenvectors with zero-structure, get the one that has
        % the largest SNR in its eigenvalue modes
        q = Q(:,I(1));
        % Enforce zero-structure
        q(M+1:end) = 0;
        % Get relative-transfer function
        h_est = q./q(param.ref);
        % Get blocking matrix
        B{n} = null(h_est(1:M)')';
        % Get fixed beamformer
        Fq{n} = h_est(1:M)/(h_est(1:M)'*h_est(1:M));
    end

    % Construct adaptive filter
    w = nan(M+L-1,N);
    for n=1:N
        w(:,n) =  pinv([B{n}*Rmm1_f(1:M,1:M,n)*B{n}' B{n}*Rmm1_f(1:M,M+1:end,n); Rmm1_f(M+1:end,1:M,n)*B{n}' Rmm1_f(M+1:end,M+1:end,n)])*[B{n}*Rmm1_f(1:M,1:M,n)*Fq{n};Rmm1_f(M+1:end,1:M,n)*Fq{n}];
    end

    % Apply filter
    % Fixed
    for n=1:N
        GEIC_GEVD.m_f(:,k,n) = Fq{n}'*signals.m_f(:,k,n) - w(:,n)'*[B{n}*signals.m_f(:,k,n);signals.l_f(:,k,n)];
        GEIC_GEVD.s_f(:,k,n) = Fq{n}'*signals.s_f(:,k,n) - w(:,n)'*[B{n}*signals.s_f(:,k,n);zeros(size(signals.l_f(:,k,n)))];
        GEIC_GEVD.n_f(:,k,n) = Fq{n}'*signals.n_f(:,k,n) - w(:,n)'*[B{n}*signals.n_f(:,k,n);zeros(size(signals.l_f(:,k,n)))];
        GEIC_GEVD.e_f(:,k,n) = Fq{n}'*signals.e_f(:,k,n) - w(:,n)'*[B{n}*signals.e_f(:,k,n);signals.l_f(:,k,n)];
    end
end

end
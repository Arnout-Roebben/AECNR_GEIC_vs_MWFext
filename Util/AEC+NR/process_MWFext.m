function MWFext = process_MWFext(signals,param)
% Processing using an extended multichannel Wiener filter (MWFext) [1,2].
%
% INPUT:
% signals      Struct     STFT-transformed input signals:
% -m_f         MXKXN      M-microphone microphone signal with K frames and N frequency bins.
% -s_f         MXKXN      M-microphone desired speech signal with K frames and N frequency bins.
% -n_f         MXKXN      M-microphone noise signal with K frames and N frequency bins.
% -e_f         MXKXN      M-microphone echo signal with K frames and N frequency bins.
% -s_f         LXKXN      M-microphone loudspeaker signal with K frames and N frequency bins.
% -rank_s      1X1        Rank to be used in the 'GEVD' procedure for the desired speech correlation matrix.
% -lambda      1X1        Exponential smoothing factor for correlaton
%                         matrix averaging, signifying the weight associated 
%                         to the previous estimate.
% -VADs        KXN        Voice activity detector (VAD). A 1 denotes
%                         desired speech activity and a 0 desired speeech inacitvity.
%
% OUTPUT
% MWFext       Struct     STFT-transformed output signals:
% -m_f         MXKXN      See INPUT
% -s_f         MXKXN      See INPUT
% -n_f         MXKXN      See INPUT
% -e_f         MXKXN      See INPUT
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

%% Initialisation
M = size(signals.m_f,1); % Number of microphones
K = size(signals.m_f,2); % Number of frames
N = size(signals.m_f,3); % Number of channels
L = size(signals.l_f,1); % Number of loudspeakers

% Preallocate memory
MWFext = struct(); % Struct to hold results
MWFext.m_f = nan(M+L,K,N); % Processed microphone signal
MWFext.s_f = zeros(M+L,K,N); % Processed desired speech signal
MWFext.n_f = zeros(M+L,K,N); % Processed near-end room directional noise signal
MWFext.e_f = zeros(M+L,K,N); % Processed echo signal

%% Processing
% Initialise correlation matrices
% Collected when desired speech is active
Rmm1_f(:,:,:,1) = repmat(1e-6*eye(M+L),[1 1 N]); 
% Collected when desired speech is inactive
Rmm0_f(:,:,:,1) = repmat(1e-6*eye(M+L),[1 1 N]); Rmm0_f(M+1:end,:,:,1) = 0; Rmm0_f(:,M+1:end,:,1) = 0;

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
                permute(zeros(size(signals.l_f(:,k,n))),[1 3 2])))*(cat(1,squeeze(signals.m_f(:,k,n)),...
                permute(zeros(size(signals.l_f(:,k,n))),[1 3 2])))');
        end
    end

    % Compute the filter
    W = nan(M+L,M+L,N);
    for n=1:N
        [~,V,Q] = updateDifferenceCorrelation(Rmm1_f(:,:,n),Rmm0_f(:,:,n),param.rank);
        l1 = diag(V'*Rmm1_f(:,:,n)*V);
        l0 = diag(V'*Rmm0_f(:,:,n)*V);
        [~,I] = sort(sum(abs(Q(M+1:end,:)),1),'ascend');
        I = sort(I(1:M),'ascend');
        Rss_est = Q(:,I(1:param.rank))*diag(max(l1(I(1:param.rank))-l0(I(1:param.rank)),0))*Q(:,I(1:param.rank))';
        W(:,:,n) = pinv(Rmm1_f(:,:,n))*Rss_est;
    end     

    % Apply the filter...
    MWFext.m_f(:,k,:) = applyFilterMultichannel(cat(1,squeeze(...
        signals.m_f(:,k,:)),permute(signals.l_f(:,k,:),[1 3 2])),W(:,:,:));
    %... to the extended desired speech signal
    MWFext.s_f(:,k,:) = applyFilterMultichannel(cat(1,squeeze(...
        signals.s_f(:,k,:)),zeros(L,N)),W(:,:,:));
    %... to the near-end room white noise signal
    MWFext.n_f(:,k,:) = applyFilterMultichannel(cat(1,squeeze(...
        signals.n_f(:,k,:)),zeros(L,N)),W(:,:,:));
    %... to the extended far-end room echo signal
    MWFext.e_f(:,k,:) = applyFilterMultichannel(cat(1,squeeze(...
        signals.e_f(:,k,:)),permute(signals.l_f(:,k,:),[1 3 2])),W(:,:,:));      
end

%% Return results
MWFext.m_f = MWFext.m_f(1:M,:,:);
MWFext.s_f = MWFext.s_f(1:M,:,:);
MWFext.n_f = MWFext.n_f(1:M,:,:);
MWFext.e_f = MWFext.e_f(1:M,:,:);

end
function VADs = VAD(x_f,sens,ref)
% Creates a binary voice activity detector (VAD) for each channel as
% abs(squeeze(s_f(ref,:,n)))> std(s_f(ref,:,n))*sens). 
%
% INPUT:
% x_f    MXKXN     M-microphone signal of K frames and N channels.
% sens   1X1       Sensitivity of the standard deviation in the VAD 
%                  formula.
% ref    1X1       Reference channel
%
% OUTPUT:
% VADs   KXN       1 denotes voice activity and 0 denotes no voice
%                  activity.
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
[~,K,N] = size(x_f); % K number of frames and N channels

%% Processing
VADs = nan(K,N); % Placeholder for VAD
for n=1:N % Loop over bins
    % Power-based VAD estimation
    VADs(:,n) = abs(squeeze(x_f(ref,:,n))) > std(x_f(ref,:,n))*sens;
end

end
function x = WOLA_synthesis(X,win,N,shift)
% Weighted overlap add (WOLA) synthesis filterbank. 
% 
% INPUT:
% X         MXKXN   Frequency matrix with K number of frames.
% win       NX1     Window.
% N         1X1     DFT size.
% shift     1X1     Frame shift.
%
% OUTPUT:
% x         TXM     Vector in time domain of length T samples.
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
M = size(X,1); % Number of microphones
% Placeholder for the output, Length of x calculated using reverse formula
% in doc STFT
x = zeros((size(X,2)-1)*shift+N,M);

%% Processing
X = cat(3,X,flip(conj(X(:,:,2:end-1)),3)); % Restore full spectrum
% Inverse discrete Fourier transform + apply window
x_full = ifft(X,N,3,'symmetric').*repmat(permute(win,[3 2 1]),...
    [M size(X,2) 1]);

% Synthesis
for l=1:size(X,2)
    x((l-1)*shift+1:(l-1)*shift+N,:) = x((l-1)*shift+1:(l-1)*shift+N,:)+...
        permute(x_full(:,l,:),[3 1 2]);
end
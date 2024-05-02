function [S, mxl] = fast_acov_spectr3(x, nfft, fs)
% [S, f] = fast_acov_spectr3(x, nfft, fs);
% Computes spectra from the FFT of cross-covariances, ignoring NaNs in the
% in the data. Returns maxlag=nfft/2 spectral estimates S.
% V.2 returns _physical_ frequencies f = fs * (1/nfft : 1/nfft : 0.5)'.
% Calls autocov function nanacov(x, maxlag),
% which ignores NaNs without interpolating across them.
% v.3 Zero pad the autocovariance before the FFT to get nfft finite spectral estimates.
%     User should zero the spectral estimates thus:
%     minimum good resolved frequency index
%     mnk = ceil(nfft/(2*mxl)); % -1+1 because truncated zero but one-based indices
%     set low-frequency unresolved spectral estimates to zero like this:
%     S(1:mnk-1) = 0; % zero all estimates before good frequencies
%
% (c) Simon de Szoeke, 2019, 2020

maxlag = nfft / 2;
ac = nanacov(x, maxlag); % no windowing necessary, returns ac(1:maxlag+1)

mxl = maxlag;
% truncate NaNs from (end of) ac before FFT; % SPdeS 2020-08-17
if any(isnan(ac))
    mxl = find(isnan(ac), 1, 'first')-1;
end

% fill low-frequency missing part of ac with zeros
Sfft = fft( [ac(1:mxl+1); zeros(nfft-2*mxl,1); ac(mxl:-1:2)], nfft);

% power spectral density in physical units
% truncate f=0 estimate, don't multiply by 2. FIXED 2020-05-20
S = 2/fs * abs( Sfft(2:maxlag+1) );

% minimum good resolved frequency index
mnk = ceil(nfft/(2*mxl)); % -1+1 because truncated zero but one-based indices
% set low-frequency unresolved spectral estimates to zero like this:
S(1:mnk-1) = 0; % zero all estimates before good frequencies

% correct but commented out to be faster:
% f = fs * (1:nfft/2)' / nfft; % frequency in physical units using sampling frequency fs, ignore f=0

% SPdeS 2020-08-17
% Note the lowest frequency can evaluate outside the function.
% Don't zero out spectral estimates from long missing lags

% % minimum resolved frequency discrete spectrum units of samples
% % set low unresolved spectral estimates to zero like this:
% S(1:mnk-1) = 0; %


% f = (1/nfft : 1/nfft : 0.5)'; % correct in _physical_ units, agrees w/ above
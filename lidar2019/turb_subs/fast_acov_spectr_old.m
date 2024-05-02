function [S, f] = fast_acov_spectr(x, nfft)
% [S, f] = fast_acov_spectr(x, nfft);
% Computes spectra from the FFT of cross-covariances, ignoring NaNs in the
% in the data. Returns maxlag=nfft/2 spectral estimates S.
% Optionally returns frequencies f = (1/nfft : 1/nfft : 0.5)'.
% Calls autocov function nanacov(x, maxlag),
% which ignores NaNs without interpolating across them.
%
% (c) Simon de Szoeke, 2019

maxlag = nfft / 2;
ac = nanacov(x, maxlag); % no windowing necessary, returns ac(1:maxlag+1)
Sfft = fft( ac([1:maxlag+1, maxlag:-1:2]), nfft);

S = abs( Sfft(2:maxlag+1) ); % ignore f=0 estimate, don't multiply by 2
f = (1:nfft/2)' / nfft; % in 1/dt sampling-frequency units, ignore f=0
%f = (1/nfft : 1/nfft : 0.5)';
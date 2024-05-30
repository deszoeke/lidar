function [P, f] = fast_acov_spectr2(x, nfft, fs)
% [S, f] = fast_acov_spectr2(x, nfft, fs);
% Computes spectra from the FFT of cross-covariances, ignoring NaNs in the
% in the data. Returns maxlag=nfft/2 spectral estimates S.
% V.2 returns _physical_ frequencies f = fs * (1/nfft : 1/nfft : 0.5)'.
% Calls autocov function nanacov(x, maxlag),
% which ignores NaNs without interpolating across them.
%
% (c) Simon de Szoeke, 2019, 2020

maxlag = nfft / 2;
ac = nanacov(x, maxlag); % no windowing necessary, returns ac(1:maxlag+1)
Sfft = fft( ac([1:maxlag+1, maxlag:-1:2]), nfft);

S = abs( Sfft(2:maxlag+1) ); % ignore f=0 estimate, don't multiply by 2. FIXED 2020-05-20
P = 2/fs * S;
% 2019-12-20 thought I should multiply by 2, but based on a wrong test in test_spectrum.m
% S = 2 * abs( Sfft(2:maxlag+1) );
f = fs * (1:nfft/2)' / nfft; % frequency in physical units using sampling frequency fs, ignore f=0
% f = (1:nfft/2)' / nfft; % in 1/dt sampling-frequency units, ignore f=0
%f = (1/nfft : 1/nfft : 0.5)'; % correct, agrees w/ above
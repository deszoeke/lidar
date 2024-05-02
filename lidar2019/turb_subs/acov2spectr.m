function [S, nu, Sfft, nufft] = acov2spectr(nfft, au)
% [S, nu] = acov2spectr(nfft, au)
% Returns the variance spectral density, the discrete Fourier transform of
% the symmetric autocovariance au(l = -maxlag:maxlag) [2*maxlag+1].
% S(nu) [nfft+1] spectral variance density
% nu[nfft+1] frequency in cycles/step [ -1/2  1/2 ]; f = nu/dt/.
% Does fortran sums so numerics don't depend on size of arrays. S can be a
% different size than au.
% For real-valued autocov, S is symmetric with useful part:
% nu(nfft/2+1:end),  2*real(S(nfft/2+1:end)
%
% Simon de Szoeke, 2018-11-29

% nfft = 128;
% N = 2*maxlag+1; % length of au
% Does fortran accumulated sums so numerics don't depend on size of arrays. S can be a
% different size than au.
maxlag = floor( length(au)/2 );
k = -(nfft/2):(nfft/2); % length nfft+1
nu = k/nfft; % frequency in cycles/step [ -1/2  1/2 ] % f = nu/dt
l = -maxlag:maxlag;
S = zeros(nfft+1,1);
for kind = 1:length(k)
    for lind = 1:length(l)
        S(kind) = S(kind) + au(lind) * exp(-2*pi*1i * nu(kind) * l(lind));
    end
end


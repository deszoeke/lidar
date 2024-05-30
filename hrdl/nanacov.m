function ac = nanacov(x, maxlag)
% ac = nanacov(x, maxlag)
% Computes the lag autocovariance ac [maxlag+1] of the vector x
% from l=0 to l=maxlag for finite data x, i.e. ignoring NaNs.
% ac(1) is the (l=0 autoco)variance.
% ac(maxlag+1) is the l=maxlag autocovariance.
%
% Performance: Comparable speed to the JIT-compiled xcov for nx=1e3, maxlag=32.
% Almost as fast as xcov (which doesn't ignore NaNs) for nx=1200, maxlag=64.
% 4x slower than xcov with nx=10^6, maxlag=32.
%
% (c) 2019, Simon de Szoeke

%isf  = isfinite(x);
%isfi = find(isf);
%xp = x - mean( x(isf) );
isfi = find(isfinite(x)); % saves memory, doesn't affect speed
xp = x - mean( x(isfi) ); % subtract one record mean
nx = length(isfi);

a = zeros(maxlag+2, 1); % last bin is a junk bin
n = zeros(maxlag+2, 1); % last bin is a junk bin
for i0 = 1:nx
    for l = 0:min(maxlag, nx-i0)
        d = min( isfi(i0+l) - isfi(i0), maxlag+1 ); % indexes one junk bin
        a(d+1) = a(d+1) + xp(isfi(i0)) * xp(isfi(i0+l));
        n(d+1) = n(d+1) + 1; % increment number of obs.
    end % l
end     % i0

% normalize and junkectomy
ac = a(1:maxlag+1) ./ n(1:maxlag+1); % biased estimator normalization
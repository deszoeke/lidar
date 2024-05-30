function [epsfilled, epx, epy] = filter_epsilon(epsilon, F ,spec, ulev, relerreps, count)
% [epsfilled, epszint, epszsmooth] = filter_epsilon(epsilon, spec, ulev, relerreps, count)
% Fill and filter TKE dissipation rate (epsilon) time series.
% Filtering expected before calculating MLD with calc_mld.m
% May not be needed for 2020 ATOMIC microDop.
% 
% Simon de Szoeke :: ATOMIC 2020 :: 2021-Apr-30

%% constants
crosscomp  = 3 / 4;
kolmogorov = 0.54; % Matches atmospheric boundary layer estimates and Sreenivasan 1995
% kolmogorov=0.5; % probably only 1 digit of precision, Sreenivasan 1995
C1prime = 4 / 3 * kolmogorov; % as in Pope eqn 6.243
factr   = crosscomp / kolmogorov; % 1/C1prime, used for my dissipation calculation

% sometimes epsilon can't be calculated properly, but it can be inferred
% to be very small, in which case
%% dissipation upper bound
% using the max spec est of the UNcompensated spectrum works well:
[smx, imx] = max(spec, [], 2);
% % epsilon defined by the 95%ile of the white noise gives too small results:
% epsnoise95 = ( factr .* (2*pi./ulev).^(2/3) .* F53(imx) .* Snoise95 ).^1.5 ;
epsUB = max(epsilon, ( factr .* (2*pi./ulev).^(2/3) .* squeeze(F(imx).^(5/3).*smx) ).^1.5  );
% Note the max of the compensated spectrum F53*spec gives too high results.
% epsilon filled with uppper bound when missing:
epsfilled = min(epsilon, epsUB); % looks good in pcolor plots
% Make missing again if rel. error is large:
epsfilled(relerreps>=1) = NaN;
epsfilled(count<700) = NaN; % too few nonnoise returns

%% interpolate dissipation gaps vertically
% epsilon==0 is a missing value, I should fix.
kmax=60; % limit MLD 
epx = epsfilled; %epsilon;
% epx(count<700 & epsilon<=0) = NaN; % filter out bad stuff at top, etc.
% epsilon = 0 is valid, i.e. too small
ii=isfinite(epx); 
for it = 1:size(epsilon,1)
    if sum(ii(it,:))>1
        epx(it,~ii(it,:))=interp1(find(ii(it,:)),epsilon(it,ii(it,:)),find(~ii(it,:)));
    end
end

%% vertically smooth
epy = runningnanmean2(epx, 3, 1, 2, 3); % runningnanmean(timeseries,window,minnobs,dim,iter)

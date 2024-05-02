function [ epsilon, sigv2, spec, maxl, Sthresh, Snoise ] = dissipation_spectr_acov3( Dv, Ulev, fs )
% dissipation_spectr_acov( Dv, fs, Urel )
% Dv   Doppler velocity (m/s)
% Urel relative wind speed for Taylor transport (m/s)
% fs   sampling frequency (s^-1)
% Compute vertical profile of dissipation from horizontal variation of
% vertical velocity using a Taylor transform at mean velocity Urel.
%
% Simon de Szoeke

% addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/stat/'));

%% turbulence constants
crosscomp  = 3 / 4;
kolmogorov = 0.54; % Matches atmospheric boundary layer estimates and Sreenivasan 1995
% kolmogorov=0.5; % probably only 1 digit of precision, Sreenivasan 1995
factr   = crosscomp / kolmogorov; % 1/C1prime, used for my dissipation calculation
% S(w; k) = C1primt * epsilon^2/3 k^-5/3
% universal constant for 2nd order structure function
% longitudinal structure function constant C2 from Pope 6.2 (p.193) after Saddoughi and Veeravalli (1994)
% C2ll   = 2.0;
% C1prime = 4 / 3 * kolmogorov; % as in Pope eqn 6.243
% factrz = 1 / C2ll;
% factrx = 3 / 4 / C2ll;

% Kolmogorov (viscous) scale for atmospheric boundary layer is
% eta=(nu^3/epsilon)^1/4; epsilon~3e-4 m^2/s^3, kinematic viscosity nu=2e-5 m^2/s
% --> eta= 2.3 mm
%     25 m = 1000 eta;  1000 m = 440 000 eta

%% Noise threshold for HildebrandSekhon3() SPdeS 2018-07-09
p1side = 0.95;
z1side = norminv( p1side );

%% w spectra window and FFT parameters used for estimating dissipation
nwin = 512; % gives 256 spectral variance esimtates, all for nonzero frequencies
lok  = 1; %  v4
hik  = 64;
nrange = 104;
kmax = 104;  %highest gate to try to calculate epsilon

F  = (1:nwin/2)'/nwin*fs; % frequencies, Hz; first is fs/nwin>0, last is Nyquist fs/2
nf = length(F);
dF = F(2) - F(1); % scalar = fs/nwin
F53   = F.^(5/3);

dtwin  = 7.5*60;         % seconds per window (7.5 min)
dswin = floor(fs*dtwin); % samples per window

% time-mean relative horizontal velocity
% Ulev = nanmean( Urel' ); % transpose

nt  = length(Dv);
if nt >= dswin / 4 % skip to next iteration if window too small
    nv               = sum( isfinite(Dv), 2 );
    [nmax, knmax]    = max( nv ); % nv is the max number of nonnoise returns, knmax is its height index
    
    % choose the highest vertical level ik that has >0.25*nmax of this coverage
    ik  = min(kmax, find(nv > 0.25*nmax, 1, 'last'));  % highest
    hk  = find(nv > 0.5*nmax, 1, 'first'); % lowest
    kk  = (max(hk-1, 1) : min(ik+1, kmax))'; % vertical indices of good data
    nk  = length(kk); % number of good levels
    w = Dv(kk,:)'; % put time dimension first
    
    % filter out extreme values
    w( abs(w-nanmean(w)) > 4*nanstd(w) ) = NaN;
        
    S = zeros( nf, 1 );
    spec = deal( NaN(nf, nk) );
    epsilon = deal( NaN(nrange,1) );
    [ Sthresh, Snoise, sigv2 ] = deal( NaN(nk, 1) );
    [knoise, maxl, countS] = deal( zeros(nk, 1) );
    
    for k = 1 : nk % loop through vertical levels
        if sum(isfinite(w(:,k))) > dswin/3 % is there enough data for spectra at this vert level
            %[S, mxl] = fast_acov_spectr3(x, nfft, fs)
            [S(:), mxl] = fast_acov_spectr3(w(:,k), nwin, fs); % Ignore Nyquist-frequency power estimate!
            % NaNs in autocov at long lags ignored. Should be more finite data 2020-07-21
            maxl(kk(k)) = mxl;  % index of longest good lag autocov
            mink = ceil(nwin/(2*mxl)); % index of first good spectral estimate
            
            if isfinite(S(1))
                spec(:, kk(k)) = S;
                [Sthresh(kk(k)), Snoise(kk(k)), knoise(kk(k))] = HildebrandSekhon3(S(mink:end), z1side);
                hfk  = min(hik, knoise(kk(k))-1); % high-frequency cutoff spectral index (inclusive)
                keep = max(lok,mink) : hfk;
                countS(kk(k)) = length(keep);
                if length(keep) < 1
                    epsilon(kk(k)) = NaN;
                else
                    vls = factr .* (2*pi / Ulev(kk(k)))^(2/3) .* ...
                        mean(F53(keep) .* (S(keep) - Snoise(kk(k)))); % simple mean
                    % mednmean(F53(keep) .* (S(keep) - Snoise(kk(k))), 5); % mean of 5 middle points
                    epsilon(kk(k)) = vls .^ 1.5; % dissipation
                end
                % calculate nonnoise variance from resolved spectrum
                sigv2(kk(k)) = dF * sum( S( 1:knoise(kk(k)) ) - Snoise( kk(k) ) );
            end % finite spectra were returned
        end % there is enough data
    end % loop vertically 

end % there is 10 min velocity data

% get rid of epsilon=0 where countS=0
epsilon(countS==0) = NaN; % shouldn't be needed, but it is



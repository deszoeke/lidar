% Calculate TKE dissipation from CSD microDop (uD) Doppler velocity spectrum
% Simon de Szoeke 2013.04.26
% revised 2018.05.29, 2018.07.25 to do HRDL data
% rev. on ATOMIC 2020 and then 2020.06.08 for ATOMIC
%
% v.7.1 Alan updated again 2021-02-24, added w CFADs, fix chunk timing glitch
% v.7 use updated netcdf files from Alan, June 2020
% v.6 fix normalization of spectrum from autocov calc.
% v.5, like v.3 but increase the window size nwin.
% v4.5 test using a window of 192 points
% 4.0 scrap gap interp; calculate spectra by FFT of autocovariance, which is insensitive to gaps.
% 3.7 compensate for spectral gain of gap interpolation
% 2.6 nwin=192, interpolate only across gaps of <=35 points


%% set environment and constants
%cd /Users/sdeszoek/Data/cruises/ATOMIC_2020/RHB/scientific_analysis/microDop
addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/stat/'));
addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/graphics/'));

rangeFuncpath = '/Users/sdeszoek/Data/cruises/ATOMIC_2020/RHB/microDop/reprocessed20200613/';
mDpath = '/Volumes/data partition on WD 25F3 SSD/Data/cruises/ATOMIC_2020/RHB/microDop/reprocessed20210224/';
% e.g. file
% EUREC4A_ATOMIC_RonBrown_lidar_microDop02_lev1_20200211_1000-20200211_1100.nc
fstem = 'EUREC4A_ATOMIC_RonBrown_lidar_microDop02_lev1_';

% quicklookpath = '/Users/sdeszoek/Data/cruises/ATOMIC_2020/RHB/microDop/quicklook';
% windprofpath = fullfile(mDpath, '/windprofiles/');
windprofpath = '/Users/sdeszoek/Data/cruises/ATOMIC_2020/RHB/microDop/windprofiles/';
samospath = '/Users/sdeszoek/Data/cruises/ATOMIC_2020/RHB/scientific_analysis/microDop';
fluxpath = '/Users/sdeszoek/Data/cruises/ATOMIC_2020/RHB/PSL/ATOMIC_2020/RHB/flux/met_sea_flux_nav/EUREC4A_ATOMIC_RonBrown_10min_nav_met_sea_flux_20200109-20200212_v1.3.nc';
P = @(x) ncread(fluxpath, x);

paren = @(x, varargin) x(varargin{:}); % for indexing or splatting arguments into functions

warning('off','MATLAB:interp1:NaNinY')

%% turbulence constants
crosscomp  = 3 / 4;
kolmogorov = 0.54; % Matches atmospheric boundary layer estimates and Sreenivasan 1995 % kolmogorov=0.5; % probably only 1 digit of precision, Sreenivasan 1995
C1prime = 4 / 3 * kolmogorov; % as in Pope eqn 6.243
factr   = crosscomp / kolmogorov; % 1/C1prime, used for my dissipation calculation
% S(w; k) = C1primt * epsilon^2/3 k^-5/3
% Universal constant for 2nd order structure function, longitudinal structure function constant C2 from Pope 6.2 (p.193) after Saddoughi and Veeravalli (1994)
C2ll   = 2.0;
factrz = 1 / C2ll;
factrx = 3 / 4 / C2ll;
% Note: Kolmogorov (viscous) scale for atmospheric boundary layer is
% eta=(nu^3/epsilon)^1/4; epsilon~3e-4 m^2/s^3, kinematic viscosity nu=2e-5 m^2/s
% --> eta= 2.3 mm
%     25 m = 1000 eta;  1000 m = 440 000 eta

%% w spectra window and FFT parameters used for dissipation
%nwin=64;
%keep=2:8; % indices of spectra to keep - not good if detrended
% v.3
% nwin = 128; % gives 64 spectral variance esimtates, all for nonzero frequencies
% nwin = 256; % gives 128 spectral variance esimtates, all for nonzero frequencies
nwin = 512; % gives 256 spectral variance esimtates, all for nonzero frequencies
% lok  = 2;
lok  = 1; % try this for v4
% hik  = 16;
hik  = 64;

fs = 2.0; %Hz
f  = (1:nwin/2)'/nwin; % discrete frequencies; variance = sum(S) * df
df = f(2) - f(1);
F  = (1:nwin/2)'/nwin*fs; % physical frequencies, Hz; first is fs/nwin>0, last is Nyquist fs/2
dF = F(2) - F(1); % scalar = fs/nwin
F53   = F.^(5/3);
dten  = 10*60; % seconds per window (10 min)
diten = floor(fs*dten); % samples per window
di    = diten;

% Noise threshold for HildebrandSekhon3() SPdeS 2018-07-09
p1side = 0.95;
z1side = norminv( p1side );

%% set up file reading
year  = 2020;
yyyy  = sprintf('%04d',year);
% files = dir( fullfile(mDpath, '/2020-*/microDop*.mat') ); % all hourly files
% files = dir( fullfile(mDpath, '*_microDop02.nc') ); % all hourly files
files = dir( fullfile(mDpath, [fstem '*.nc']) ); % all hourly files
dnum_file = zeros( length(files), 1 );
for i = 1:length( files ) % 881 files
    dnum_file(i) = datenum( files(i).name(length(fstem)+(1:13)), 'yyyymmdd_HHMM' );
end
istart = find( dnum_file >= datenum( year, 01, 07 ), 1, 'first');
istop  = length( files );

% sugar, low winds, dust case:
% Jan 29 - Feb 3
% istart = find( dnum_file >= datenum( year, 01, 29 ), 1, 'first');
% istop  = min( length( files ), find( dnum_file < datenum( year, 02, 04 ), 1, 'last') );
% note in port Jan 26-28 2020.

%% load height coordinate from Doppler file
retrfile = fullfile(files(istart).folder, files(istart).name);
% load( retrfile ); % --> mdData
% zData = mdData.data.rangeArr(1, :); % zData = range variable for Doppler w
D = @(x) ncread( retrfile, x );
zData = D('range') + D('scannerHeight'); % SPD added scannerHeight 2020-07-02
iz = isfinite(zData);
zData = zData(iz); % truncate missings (NaN) at end

%% load range function for SNR to range-correct backscatter intensity
load(fullfile(rangeFuncpath, 'rangeFunc_ATOMIC.mat')); % rangeFun

%% load ship-relative wind profiles
% need to update wind profiles
% wnd = load(fullfile( windprofpath, 'windprofs.mat' )); % microDop ship-relative wind profiles
% wnd. twind wind uwind vwind speedwind --> uvad, vvad
% new NetCDF wind profiles from Alan, 2020 Jun 10:
windfile = fullfile(windprofpath, 'ATOMIC_RHB_mDop_windProf.nc');
wnd.zwind = ncread(windfile, 'height');
height_mid_wind = 0.5*( wnd.zwind(1:end-1) + wnd.zwind(2:end) );
dz_wind = diff( wnd.zwind );
dz_smth = dz_wind;
% dz_smth(2:end) =  (height_mid_wind(2:end)/2.2e-3).^(1/4.2);
wnd.timewind = ncread(windfile, 'yDay'); % decimal day of year
wnd.speed = ncread(windfile, 'speed');
wnd.direction = ncread(windfile, 'direction');
wnd.uwind = -wnd.speed .* sind(wnd.direction);
wnd.vwind = -wnd.speed .* cosd(wnd.direction);
shear2 = (diff(wnd.uwind, 1, 1)./dz_smth).^2 + (diff(wnd.vwind, 1, 1)./dz_smth).^2;
% 2020-06-10: Winds are earth-relative wind profiles now!
% ship maneuvers not removed!!

%% load ship horizontal velocity --need to update SAMOS winds for leg2
load(fullfile( samospath, 'samos.mat' ))  % samos.md0, samos.ushp, samos.vshp
%% Relative winds
wnd.urel = wnd.uwind - ...
        interp1(samos.md0(:), samos.ushp, datenum(2020, 0, wnd.timewind))';
wnd.vrel = wnd.vwind - ...
        interp1(samos.md0(:), samos.vshp, datenum(2020, 0, wnd.timewind))';

%{
b2rbrewcmap(8)
subplot(2,1,1)
pcolor(wnd.timewind, wnd.zwind/1e3, wnd.uwind); shading flat
caxis([-16 16]); ylim([0 2.5])
subplot(2,1,2)
pcolor(wnd.timewind, wnd.zwind/1e3, wnd.vwind); shading flat
caxis([-16 16]); ylim([0 2.5])
subplot(2,1,2)
pcolor(wnd.timewind, height_mid_wind/1e3, log10(shear2)); shading flat
caxis([-8 -2.5]); ylim([0 2.5])
%}

% interpolate relative wind to w-levels
isf = isfinite(wnd.zwind);
[uwind_w, vwind_w] = deal(NaN(size(wnd.uwind, 2), length(zData)));
for i = 1:size(wnd.uwind, 2) % for all times
    isfz = isf & isfinite(wnd.urel(:,i));
    if any(isfz)
        uwind_w(i,:) = interp1(wnd.zwind(isfz), wnd.urel(isfz,i), zData(:));
        vwind_w(i,:) = interp1(wnd.zwind(isfz), wnd.vrel(isfz,i), zData(:));
    end
end
isft = isfinite(wnd.timewind);
intu = @(t) interp1(wnd.timewind(isft), uwind_w(isft,:), t);
intv = @(t) interp1(wnd.timewind(isft), vwind_w(isft,:), t);
% intu(ydmid-datenum(year,0,0))
% intv(ydmid-datenum(year,0,0))

%% turbulence analysis
% allocate output data variables
nhr = 24 * 6;
n10min = nhr * 6;
kmax = 60; % find(zData<2e3,1,'last')
ndk = 5;
S = NaN(size(F));
count   = zeros(n10min,kmax);
epsilon = NaN(n10min,kmax);
ulev    = NaN(n10min,kmax);
% epsx=NaN(n10min,kmax);
% epsz=NaN(n10min,ndk);
ydayw10 = NaN(n10min,1);
sigv2   = NaN(n10min,kmax);
spec    = NaN(n10min,floor(nwin/2),kmax);
knoise  = zeros(n10min, kmax, 'uint8');
Snoise  = zeros(n10min, kmax);
Sthresh = zeros(n10min, kmax);

% bin edges for velocity CFADs
wedges = [-Inf -7:0.1:7 Inf];
% velocity histogram for each 10-min chunk, to be filled by cfad(w, edges)
Hw = zeros( n10min, kmax, length(wedges)-1 );

%% loop and compute dissipation, or load it from a prev file
recompute_dissipation = false;% || true;
if recompute_dissipation
    %% loop over hourly files
    % handle vertical stares split across hourly files
    itmax   = 4000;
    itrmax  = 12700;
    yday    = NaN(itmax, 1);
    W       = NaN(itmax, kmax);
    
    bigind  = 0; % counter for filling output arrays
    % must zero bigind before starting ifile loop!
    
    for ifile = istart : istop
        fprintf(1,'%s\n', files(ifile).name) % microDop file
        retrfile = fullfile(files(ifile).folder, files(ifile).name);
        %         load( retrfile ); % --> mdData
        %         yDay = mdData.decTime; % matlab datenum
        D = @(x) ncread( retrfile, x ); % data loader
        % D('time') is now seconds of year
        yDay = datenum( D('year'), 0, 1, 0, 0, D('time') ); % matlab datenum
        % Jan 1 error fixed in 7.1! SPd 2021 Apr 30
        
        % range function for each leg
        if yDay(1) < rangeFun(2).startTime % matlab datenums
            rngFunc = rangeFun(1).blendRngFnc; % linear scale
        else
            rngFunc = rangeFun(2).blendRngFnc; % linear scale
        end
        
        % get motion variables
        % inuVelDown  = paren( D('inuVelDown'), ':'); % already synchronized
        %     VelCor = mdData.header.inu2ScanOffset.applyVelCor; % Alan?
        %     inuCompTime = mdData.header.mcData.inuCompTime(:);
        %     inuUtcTime  = mdData.header.mcData.inuUtcTime(:);
        %     inuGpsTime  = mdData.header.mcData.inuGpsTime(:);
        %     inuVelBodyZ = mdData.header.mcData.inuVelBodyZ(:);
        % velocity data
        %         D = mdData.data;
        % housekeeping, timekeeping, identify scan
        %         sz = size(D.decTimeArr);
        mcVelData = D( 'mcVelData' )'; % already corrected
        sz = size( mcVelData );
        % dvec = datevec( yDay );
%         scanminute = dvec(:,5) < 1; % 1 minute PPI at start of each hour
%         nt = sum( ~scanminute );
        %         dvec = reshape( datevec(D.decTimeArr), [sz, 6]);
        %         scanminute = dvec(:,:,5) < 1; % :00 - :01 is a 65-degree-elev. azimuthal scan
        
        % replace above filtering condition with scanType==0 designation for staring
        scanminute = D('scanType')~=0; % 0=vertical stare, 7=PPI
        % also blank out the discontiguous bit of data before the PPI
        scanminute(1:find(scanminute, 1, 'first')) = true;
        nt = sum( ~scanminute );
        
        % quantity proportional to backscatter --not needed for noise filtering
        wbSnr = paren(D('snrData'), iz, ':')'; % --> dB (time, height)
        % bscat2 assumes r^-2 range function; SPD used before Alan's range function.
        bscat2 = 10.^(wbSnr/10) .* zData'.^2; % linear units
        % Alan's "overlap" range function, 2020-07-03
        bscat = 10.^(wbSnr/10) ./ rngFunc(1:210); % linear units
        % Alan makes range corrected intensity (RCI):
        % rci(jk,:) = 10*log10(accum.data_wbSnr(dex,jk)/rngFunc (jk))';
        % nt = length(yDay); % will truncate to exclude VAD scan

        %% condition data on good returns
        % exclude noise, too high, azimuthal scans
        % depends on SNR, not range corrected intensity!
        % SPD eyeballed threshold during ATOMIC.
        zap = 10.^(wbSnr/10) < 3e-3 & D('recordType'); %recordType==0 -> noise
        % mask selects valid returns
        mask = (~zap  &  zData' < 2e3) & ~scanminute; %
                                         % staring after scan: scanminute defined by scanType~=0
        % repmat(~scanminute, [1 sz(2)]);
        % estimate valid velocity range
        vrng = min(5, max(abs(quantile( mcVelData(mask), [0.08 0.92] ))));
        nanzap = zeros( size(mask) );
        nanzap( ~mask ) = NaN;
        %bscat = bscat + nanzap;
        
        % motion-adjusted Doppler velocity, shift INU by 1 time step,
        % cut off scanminute
        %         platvel = inuVelDown([1 1:end-1]);
        % platvel = inuVelDown;
        W =      nanzap(~scanminute,1:kmax) ...
            + mcVelData(~scanminute,1:kmax) ; %...
        % platform velocity is already added in mcVel
%             -  repmat( platvel(~scanminute), [1, kmax] );
        
        % NOT USED: Throw out extreme velocities based on standard deviation
        % throww = abs(W) > 4.0*nanstd(W(:)); % not thrown out yet!
        % W(abs(W) > 4.0*nanstd(W(:))) = NaN;
        
        %% find chunk start and end indices
        % Locked to 10-min intervals and to file boundaries.
        % Not using rolling (e.g. FIFO) registers to handle 10-min chunks across files.
        % from diff(yday) gaps of more than 5 minutes
        %         dift    = find( diff(yDay( ~scanminute )) > 5/(24*60) );
        
        % this doesn't work at 00:59:59, for example!...
        %  minute = dvec( ~scanminute, 5 );
        %  second = dvec( ~scanminute, 6 );
        %  chunkind = floor( ( minute + second/60 ) / 10 ); % 0-5 index of each chunk
        % Fix timing bug SPd 2021 Apr 30
        chunkind = floor( ( yDay-floor(yDay(1)) ) * 144 ); % 0-5 index of each 10 min chunk
        % for calculating end chunks that cross the hour
        % chunkind( find(chunkind>=5, 1, 'last') : end ) = 5;
        uch = sort( unique( chunkind ) );
        nchunks = length( uch );
        
        %% loop through chunks and do analysis
        for ind = 1 : nchunks % index of the 10 minute interval in the hour
            chii = chunkind == uch(ind);
            i10  = find( chii );
            if ~isempty(i10)
                starts(ind) = i10( 1   );
                ends(ind)   = i10( end );
                nt  = length( i10 );
                %if nt >= di / 4 % skip to next iteration if window too small
                
                % choose the highest vertical level ik that has >0.25*nmax of this coverage
                nw              = sum( isfinite(W(i10,:)) );
                [nmax, knmax]   = max( nw ); % nw is the max number of nonnoise returns, knmax is its height index
                ik  = min(kmax, find(nw>0.25*nmax, 1, 'last'));  % highest
                hk  = find(nw>0.5*nmax, 1, 'first'); % lowest
                kk  = (max(hk-1, 1) : min(ik+1, kmax))';
                nk  = length(kk);

                if ( nt >= di / 4 ) && ( nk > 0 ) % skip to next iteration if window too small or no returns                    
                    bigind          = bigind + 1;          % increment index
                    ydstart(bigind) = yDay( starts(ind) );
                    ydend(bigind)   = yDay(   ends(ind) );
                    ydayw10(bigind) = yDay( starts(ind) ); % current time of this 10-minute window
                    
                    % subset
                    w10 = W(i10, kk);
                    
                    % CFAD w histograms; Hw( bigind, height, bin )
                    Hw( bigind, kk, : ) = cfad( w10 , wedges );
                    
                    % filter out extreme values -- already done on hourly nanstd
                    w10( abs( w10 - nanmean(w10) ) > ( 4 * nanstd(w10) ) ) = NaN;
                    iwt = isfinite( w10 );
                    count( bigind, kk ) = sum( iwt );
                    
                    % interpolate horizontal wind to the middle of each chunk
                    ydmid = nanmedian(yDay(i10)) ;
%                     iiw   = isfinite(wnd.timewind);
%                     jjw   = isfinite(wnd.zwind   );
%                     % wind from microDopp to w levels and chunk times
%                     ulevi = interp2(wnd.timewind(iiw), wnd.zwind(jjw), wnd.uwind(jjw,iiw), ydmid-datenum(year, 0, 0), zData(1:kmax), 'linear');
%                     vlevi = interp2(wnd.timewind(iiw), wnd.zwind(jjw), wnd.vwind(jjw,iiw), ydmid-datenum(year, 0, 0), zData(1:kmax), 'linear');
%                     % ship relative wind components
%                     ulevr = ulevi - interp1( samos.md0, samos.ushp, ydmid );
%                     vlevr = vlevi - interp1( samos.md0, samos.vshp, ydmid );
                    ulevr = intu(ydmid-datenum(year,0,0));
                    vlevr = intv(ydmid-datenum(year,0,0));
                    % ship-relative speed for Taylor transform
                    Ulev  = sqrt( ulevr.*ulevr + vlevr.*vlevr );
                    
                    for k = 1 : nk % loop through vertical levels
                        if count(bigind,kk(k)) > di / 4
                            % calculate the spectrum
                            S(:) = fast_acov_spectr2(w10(:,k), nwin, fs); % Ignore Nyquist-frequency power estimate!
                            if isfinite(S(1))
                                spec(bigind,:,kk(k)) = S;
                                % estimate dissipation from the spectrum
                                [Sthresh(bigind, kk(k)), Snoise(bigind, kk(k)), knoise(bigind, kk(k))] = HildebrandSekhon3(S, z1side);
                                hfk  = min(hik, knoise(bigind, kk(k))-1); % high-frequency cutoff spectral index (inclusive)
                                keep = lok:hfk;
                                if length(keep) < 1
                                    epsilon(bigind,kk(k)) = NaN;
                                else
                                    vls = factr .* (2*pi/Ulev(k))^(2/3) .* ...
                                        mednmean(F53(keep) .* (S(keep) - Snoise(bigind,kk(k))), 5); % mean of 5 middle points
                                    epsilon(bigind,kk(k)) = vls .^ 1.5; % dissipation
                                end
                                % calculate nonnoise variance from resolved spectrum
                                sigv2(bigind,kk(k)) = df * sum( S(1:knoise(bigind,kk(k))) - Snoise(bigind,kk(k)) ); % dF->df fixes discrete spectrum interpretation
                                ulev (bigind,kk(k)) = Ulev(k);
                            end % finite spectra were returned
                        end % there is data for spectra
                    end % k level where there is data
                    
                end % there is 10 min Doppler velocity data
            end % there is data in a chunk
        end % 10 min chunks
    end     % microDop file
    
    %% truncate output data
    ydayw10 = ydayw10(1:bigind);
    ydstart = ydstart(1:bigind);
    ydend   = ydend  (1:bigind);
    epsilon = real(epsilon); % makes imaginary epsilon (from S<Snoise) into zeros
    epsilon = epsilon(1:bigind,:);
    sigv2   = sigv2  (1:bigind,:);
    count   = count  (1:bigind,:);
    spec    = spec   (1:bigind,:,:);
    knoise  = knoise (1:bigind,:);
    Snoise  = Snoise (1:bigind,:);
    Sthresh = Sthresh(1:bigind,:);
    ulev    = ulev   (1:bigind,:);
    Hw      = Hw     (1:bigind,:,:);
    
    % epsilon = 0 for most significant spectral estimates that don't exceed noise.
    
    %% nondimensional spectra scaled by Kolmogorov power spectral density (dissipation*kinvisc^5)^1/4
    % for cloud top
    kinvisc = 1.63e-5; %m^2/s, at ~900 hPa
    Kpsd = sqrt(sqrt( epsilon *  (kinvisc^5) ));  % Kolmogorov power spectral density (dissipation*kinvisc^5)^1/4  [U^2 / L^-1]
    Kks  = sqrt(sqrt( epsilon ./ (kinvisc^3) ));  % Kolmogorov wavenumber (dissipation*kinvisc^3)^1/4  [L^-1]
    eta = 1 ./ Kks;                            % Kolmogorov length scale
    ulevx = reshape(ulev,[size(ulev,1) 1 size(ulev,2)]);
    Fx = reshape(repmat(F,[1 size(ulevx,1)])',[size(ulevx,1),length(F),1]);
    wavenumber = 2 * pi * Fx ./ ulevx;
    kscaled  = wavenumber ./ reshape(Kks,[size(ulev,1) 1 size(ulev,2)]); % nondimensional wavenumber is scaled by Kolmogorov wavenumber
    speck    = spec     .* ulevx/(2*pi);   % wavenumber spectrum (spec is frequency spectrum!)
    Sfscaled = spec     ./ reshape(Kpsd,[size(ulev,1) 1 size(ulev,2)]) ;  % nondimensional (TKE) power frequency-spectral density scaled by Kolm'v PSD
    Skscaled = Sfscaled .* ulevx/(2*pi);  % nondimensional (TKE) power wavenumber-spectral density scaled by Kolm'v PSD
    
    %% cut off noise
    Scut = Skscaled;
    kx   = max(1,knoise);
    for     i = 1 : size(knoise,1)
        for j = 1 : size(knoise,2)
            Scut(i, double(kx(i,j)):end, j) = NaN;
        end
    end
    Snoisex = reshape(Snoise,[size(ulev,1) 1 size(ulev,2)]);
    Scutadj = Scut     - Snoisex .* ulevx/(2*pi); % subtract noise
    Sadj    = Skscaled - Snoisex .* ulevx/(2*pi);
    
    %% compute variance of compensated spectrum on every level
    Scomp = Scut(:,:,:).*kscaled(:,:,:).^(5/3);
    nt = size(Skscaled,1); nspec = nwin/2;
    Svar = NaN(kmax,1);
    Nvar = zeros(kmax,1);
    for lev = 1:kmax
        Nvar(lev) = sum( paren( isfinite( Scomp(:,:,lev) ), ':' ) );
        Svar(lev) = nanvar( paren( Scomp(:,:,lev), ':' ) );
    end
    
    %% error calculations bootstrapped from variance of compensated spectrum
    % compute a pooled variance:
    % then Sum(i=1:k)( (n_i-1)*s_i^2 ) / Sum(n_i) - k
    % k is the number of dissipation estimates
    % first compute variance of spectral estimates over each spectrum
    s2 = nanvar( Scomp, 0, 2 );
    nspecvar = sum( isfinite(Scomp), 2 );
    ps2 = nansum( max(0, nspecvar - 1).*s2, 1 ) ./ nansum( max(0, nspecvar - 1), 1 ); % pooled variance
    % Pooled variance has fewer degrees of freedom overall bc. it take into
    % account that the mean dissipation is subtracted from each compensated
    % spectrum.
    sm = nanmean( reshape( Scomp, [size(Scomp,1)*size(Scomp,2), size(Scomp,3)] ) );
    
    % loglog(Svar, zData(1:kmax))
    % hold on
    % loglog(ps2(:), zData(1:kmax))
    % plot(0.25*([400 2000]/400).^-(.27),[400 2000])
    % Variance decreases with height. Pooled variance is less sensitive to height.
    % Why does variance _decrease_ with height?
    % Number of observations decreases with height. Does this decrease the pooled
    % variance among the spectral estimates?
    % ps2 is an estimate of the population variance for each spectral estimate
    % of the dissipation. The variance of N samples is ps2/N
    % Conservatively use as the nominal dissipation of a spectral estimate:
    ps2bar = 0.3;
    sbar = nanmedian(sm); % 0.78 nominal
    % relative error if compensated spectral estimates
    sigsos = sqrt(ps2(:)) ./ sm(:); % std(Scomp) / mean(Scomp) = 0.68
    sigepsoeps = sigsos.^(3/2); % ~0.55
    relerreps  = reshape(sigepsoeps,[1 1 kmax]) ./ nspecvar.^(3/4);
    % the mean of the middle 5 points is less variable than the mean
    
    %% save data
    % save uD_turb7p1_512n.mat ydayw10 ydstart ydend F spec epsilon sigv2 count ulev knoise Snoise Sthresh Scomp relerreps nspecvar zData wedges Hw
else
    load uD_turb7p1_512n.mat
end

% redo Kolmogorov spectrum calcs after saving --a bit slow
%% nondimensional spectra scaled by Kolmogorov power spectral density (dissipation*kinvisc^5)^1/4
% for cloud top
kinvisc = 1.63e-5; %m^2/s, at ~900 hPa
Kpsd = sqrt(sqrt( epsilon *  (kinvisc^5) ));  % Kolmogorov power spectral density (dissipation*kinvisc^5)^1/4  [U^2 / L^-1]
Kks  = sqrt(sqrt( epsilon ./ (kinvisc^3) ));  % Kolmogorov wavenumber (dissipation*kinvisc^3)^1/4  [L^-1]
eta = 1 ./ Kks;                            % Kolmogorov length scale
ulevx = reshape(ulev,[size(ulev,1) 1 size(ulev,2)]);
Fx = reshape(repmat(F,[1 size(ulevx,1)])',[size(ulevx,1),length(F),1]);
wavenumber = 2 * pi * Fx ./ ulevx;
kscaled  = wavenumber ./ reshape(Kks,[size(ulev,1) 1 size(ulev,2)]); % nondimensional wavenumber is scaled by Kolmogorov wavenumber
speck    = spec     .* ulevx/(2*pi);   % wavenumber spectrum (spec is frequency spectrum!)
Sfscaled = spec     ./ reshape(Kpsd,[size(ulev,1) 1 size(ulev,2)]) ;  % nondimensional (TKE) power frequency-spectral density scaled by Kolm'v PSD
Skscaled = Sfscaled .* ulevx/(2*pi);  % nondimensional (TKE) power wavenumber-spectral density scaled by Kolm'v PSD

%% cut off noise
Scut = Skscaled;
kx   = max(1,knoise);
for     i = 1 : size(knoise,1)
    for j = 1 : size(knoise,2)
        Scut(i, double(kx(i,j)):end, j) = NaN;
    end
end
Snoisex = reshape(Snoise,[size(ulev,1) 1 size(ulev,2)]);
Scutadj = Scut     - Snoisex .* ulevx/(2*pi); % subtract noise
Sadj    = Skscaled - Snoisex .* ulevx/(2*pi);

%% compute variance of compensated spectrum on every level
Scomp = Scut(:,:,:).*kscaled(:,:,:).^(5/3);
nt = size(Skscaled,1); nspec = nwin/2;
Svar = NaN(kmax,1);
Nvar = zeros(kmax,1);
for lev = 1:kmax
    Nvar(lev) = sum( paren( isfinite( Scomp(:,:,lev) ), ':' ) );
    Svar(lev) = nanvar( paren( Scomp(:,:,lev), ':' ) );
end

%% error calculations bootstrapped from variance of compensated spectrum
% compute a pooled variance:
% then Sum(i=1:k)( (n_i-1)*s_i^2 ) / Sum(n_i) - k
% k is the number of dissipation estimates
% first compute variance of spectral estimates over each spectrum
s2 = nanvar( Scomp, 0, 2 );
nspecvar = sum( isfinite(Scomp), 2 );
ps2 = nansum( max(0, nspecvar - 1).*s2, 1 ) ./ nansum( max(0, nspecvar - 1), 1 ); % pooled variance
% Pooled variance has fewer degrees of freedom overall bc. it take into
% account that the mean dissipation is subtracted from each compensated
% spectrum.
sm = nanmean( reshape( Scomp, [size(Scomp,1)*size(Scomp,2), size(Scomp,3)] ) );

% loglog(Svar, zData(1:kmax))
% hold on
% loglog(ps2(:), zData(1:kmax))
% plot(0.25*([400 2000]/400).^-(.27),[400 2000])
% Variance decreases with height. Pooled variance is less sensitive to height.
% Why does variance _decrease_ with height?
% Number of observations decreases with height. Does this decrease the pooled
% variance among the spectral estimates?
% ps2 is an estimate of the population variance for each spectral estimate
% of the dissipation. The variance of N samples is ps2/N
% Conservatively use as the nominal dissipation of a spectral estimate:
ps2bar = 0.3;
sbar = nanmedian(sm); % 0.78 nominal
% relative error if compensated spectral estimates
sigsos = sqrt(ps2(:)) ./ sm(:); % std(Scomp) / mean(Scomp) = 0.68
sigepsoeps = sigsos.^(3/2); % ~0.55
relerreps  = reshape(sigepsoeps,[1 1 kmax]) ./ nspecvar.^(3/4);
% the mean of the middle 5 points is less variable than the mean

% the order of the times is getting screwed up each hour - fixed in v.7.1
% [yday10, ord10] = sort( ydayw10 );  % re-sort times

% compute vertically filled, interpolated, and smoothed epsilon
[epsfilled, epszint, epszsmooth] = filter_epsilon(epsilon, F, spec, ulev, relerreps, count);
% mixed layer depth from smoothed epsilon
[ mld, kmld0, kmld1, epsbar, epstop, mldgrad, kgrad ] = calc_mld( epsilon, epszsmooth, zData, 3.0 );

% w histogram example
it = 1000; iz = 30;
bar(wedges(2:end-2), squeeze( Hw(it, iz,2:end-1) ), 1, 'linestyle','none')
xlim([-3 3]); xlabel('vertical velocity');
set(gca,'fontsize',14)
title([datestr(ydayw10(it),'yyyy-mm-dd HH:MM') sprintf('   %.0f', zData(iz)) ' m'],...
        'fontweight','norm')
    
% spectrogram
lev = 18; % 633 m, near cloud base
b2rbrewcmap(12)
pcolor(ydayw10(1000:2000), log10(F), log10(squeeze(spec(1000:2000,:,lev)))'); shading flat; datetick('x')
% pcolor(ydayw10, log10(F), log10(squeeze(spec(:,:,lev)))'); shading flat; datetick('x')
xlabel('2020 date'); ylabel('log10(f)')

% plot epsilon
clf
subplot(2,1,1);
% zap = zeros(size(epsilon));
% zap (epsilon == 0) = NaN;
pcolor(ydayw10, zData(1:kmax)/1e3, log10(epszint(:,:)')); shading flat;
datetick('x');
cb = colorbar;
set(gca,'color',0.7+[0 0 0], 'tickdir','out', 'fontsize', 14); caxis([-6 -3])
axis tight
ylabel('height (km)'); xlabel('2020 date'); set(gca,'xminortick','on')
% saveas(gcf,'epsilon_atomic.eps','epsc')
% saveas(gcf,'epsilon_atomic.png','png')

% zoom in on DWL days
xlim([datenum(2020,01,30), datenum(2020,02,05)])
ylim([0 1.5])
datetick('x','mm/dd','keeplimits');

% zoom in on sugar clouds and weak wind at beginning of leg 2
xlim([datenum(2020,01,29), datenum(2020,02,03)])
ylim([0 2]); yticks(0:0.5:2)
datetick('x','mm/dd','keeplimits');

% overplot wind speed and buoyancy flux
hold on
plot(datenum(2020,0,0,0,0,P('time')), P('wspd_sfc')/10, 'k', 'linewidth',1.4)
[bflx, bflx_s] = buoyancy_flux(-P('hs_bulk'), -P('hl_bulk'), P('tair'), P('qair'), P('psealevel'));
plot(datenum(2020,0,0,0,0,P('time')), bflx*1e3, 'b', 'linewidth',1.4)
title({'low-wind, sugar, dust period'
       'wind speed/(10 m s^{-1}) (black), buoyancy flux/(10^3 m^2 s^{-3}) (blue)'},...
       'fontweight','norm')
ax = gca; ax.Position = [1 0 1 1] .* ax.Position + [0 0.4 0 0];
cb.Label.String='log_{10} TKE dissipation'; cb.Label.FontSize=14;
% saveas(gcf,'epsilon_atomic_sugar.eps','epsc')
% saveas(gcf,'epsilon_atomic_sugar.png','png')


%% plot scaled not by energy containing scale or by noise, but by inferred
% Kolmogorov scales which depend only on dissipation and viscosity
lev=18;
clf
% loglog(permute(kscaled,[2 1 3]),permute(Skscaled,[2 1 3]),'.','markersize',1,'color',0.7*[1 1 1])
loglog(kscaled(:,:,lev)',Skscaled(:,:,lev)','.','markersize',0.5,'color',0.7*[1 1 1]) % all spectral estimates
hold on
loglog(kscaled(:,lok:end,lev)',Scut(:,lok:end,lev)','b.','markersize',1)
% loglog(kscaled(:,4:end)',Scutadj(:,2:end)','c.','markersize',1) % plots all levels
plot([1e-6 1e-2],C1prime*[1e-6 1e-2].^(-5/3),'r-','linewidth',1)
set(gca,'fontsize',16,'xlim',[1e-6-10*eps 1e-1])
title('ATOMIC 2020, 634 m Kolmogorov-scaled w wavenumber spectra', 'fontweight','normal')
xlim([1e-6 1e-2])
ylabel('power spectral density   S_{ww}(k)/(\epsilon\nu^5)^{1/4}')
xlabel('wavenumber   k\eta')

% plot compensated spectrum
% clf
loglog(kscaled(:,:,lev)',Skscaled(:,:,lev)'.*kscaled(:,:,lev).^(5/3)','.','markersize',1,'color',0.7*[1 1 1])
hold on
loglog(kscaled(:,2:end,lev)',Scut(:,2:end,lev)'.*kscaled(:,2:end,lev).^(5/3)','b.','markersize',1)
plot([1e-6 1e-2],C1prime*[1 1],'r-')
saveas(gcf, 'Kolmogorov_spectrum.png')

%% plot dissipation time-height profile
clf
subplot(2,1,1)
pcolor(ydayw10, (zData(1:kmax)-2.4)/1e3, log10(epsilon(:,:))'); shading flat; caxis([-4.5 -2.5]); colorbar;
datetick('x')
xlim(datenum(2020, 01, [08 15])); ylim([0 1.2])
caxis([-4.5 -2.5])
ylabel('height (km)'); xlabel(yyyy); set(gca,'fontsize',14, 'color',0.7+[0 0 0])
set(gca,'fontsize',16, 'tickdir','out')
title('TKE dissipation (m^2 s^{-3})', 'fontweight','normal')
saveas(gcf, 'atomic_dissipation.png')

%% plot the number of ?spectral? estimates used in each dissipation estimate
plot((1:sum(nspecvar(:)>0))/sum(nspecvar(:)>0), sort(nspecvar(nspecvar>0)))
hold on
% expand 0.990 1.000 over 0-1
plot((1:sum(nspecvar(:)>0))/sum(nspecvar(:)>0)*100 -99, sort(nspecvar(nspecvar>0)))
xlim([0 1])
ylabel('dof'); xlabel('fraction of dissipation estimates'); set(gca,'fontsize', 14)
% print -depsc dof_dissipation.eps



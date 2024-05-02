% Calculate TKE dissipation from HRDL Doppler velocity spectrum
% Simon de Szoeke 2013.04.26
% revised 2018.05.29, 2018.07.25 to do HRDL data
%
% v.5, like v.3 but increase the window size nwin.
% v4.5 test using a window of 192 points
% 2.6 nwin=192, interpolate only across gaps of <=35 points
% 3.7 compensate for spectral gain of gap interpolation
% 4.0 scrap gap interp; calculate spectra by FFT of autocovariance, which is insensitive to gaps.

%% set environment and constants
%cd /Users/sdeszoek/Data/cruises/DYNAMO_2011/Revelle/HRDL
addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/stat/'));
addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/graphics/'));

warning('off','MATLAB:interp1:NaNinY')
machine=char(java.net.InetAddress.getLocalHost.getHostName); % doesn't work with floating unnamed IP nodenames!
if regexp(machine,'squall')==1
    stem = '/Users/sdeszoek/Data/cruises/DYNAMO_2011/Revelle/HRDL';
    way_DopplerData = '/Users/sdeszoek/Data/cruises/DYNAMO_2011/Revelle/HRDL/DopplerData/bscan_netCdf'; % all HRDL data
    way_RadarData   = '/Users/sdeszoek/Data/cruises/DYNAMO_2011/Revelle/wband/Raw/all';
else
    stem = '/Users/sdeszoek/Data/cruises/DYNAMO_2011/Revelle/HRDL';
    way_DopplerData = '/Users/sdeszoek/Data/cruises/DYNAMO_2011/Revelle/HRDL/DopplerData/bscan_netCdf';
    way_RadarData   = '/Users/sdeszoek/Data/cruises/DYNAMO_2011/Revelle/radar/wband/Raw/all';
end

year  = 2011;
yyyy  = sprintf('%04d',year);
files = dir(fullfile(way_DopplerData,'*_hrdl_zen.nc'));
istart = 1;
while datenum(files(istart).name(1:11),'yyyymmdd_hh') < datenum(2011,0,274) %281
    istart = istart + 1; % I tried, but can't vectorize this!
end

retrfile = fullfile(way_DopplerData,files(istart).name);
zData    = mean(ncread(retrfile,'zData'),2); % zData constant not same as horizontal wind height
radarfiles = dir(fullfile(way_RadarData,'20112810000MMCRMom.nc'));
kradarmax = 66;
hradar = double( ncread( fullfile(way_RadarData, radarfiles(1).name), 'Heights', [1 1], [kradarmax 1] ) );
noisefloor = max(-45,-35.5+20*log10(hradar)); % clouds have reflectivity factor Z above this level

% turbulece constants
crosscomp  = 3 / 4;
kolmogorov = 0.54; % Matches atmospheric boundary layer estimates and Sreenivasan 1995
% kolmogorov=0.5; % probably only 1 digit of precision, Sreenivasan 1995
C1prime = 4 / 3 * kolmogorov; % as in Pope eqn 6.243
factr   = crosscomp / kolmogorov; % 1/C1prime, used for my dissipation calculation
% S(w; k) = C1primt * epsilon^2/3 k^-5/3
% universal constant for 2nd order structure function
% longitudinal structure function constant C2 from Pope 6.2 (p.193) after Saddoughi and Veeravalli (1994)
C2ll   = 2.0;
factrz = 1 / C2ll;
factrx = 3 / 4 / C2ll;

% Kolmogorov (viscous) scale for atmospheric boundary layer is
% eta=(nu^3/epsilon)^1/4; epsilon~3e-4 m^2/s^3, kinematic viscosity nu=2e-5 m^2/s
% --> eta= 2.3 mm
%     25 m = 1000 eta;  1000 m = 440 000 eta

% Noise threshold for HildebrandSekhon3() SPdeS 2018-07-09
p1side = 0.95;
z1side = norminv( p1side );

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
F  = (1:nwin/2)'/nwin*fs; % frequencies, Hz; first is fs/nwin>0, last is Nyquist fs/2
dF = F(2) - F(1); % scalar = fs/nwin
F53   = F.^(5/3);
dten  = 10*60; % seconds per window (10 min)
diten = floor(fs*dten); % samples per window
di    = diten;

%% load wind profiles
%if regexp(machine,'squall')==1
ncf = fullfile(stem,'windProf.nc');
yDay_wind   = ncread(ncf,'yDay');
height_wind = ncread(ncf,'height');
speed       = ncread(ncf,'speed');
direction   = ncread(ncf,'direction');
%elseif regexp(machine,'fog')==1
%end
uvad = -speed .* sind(direction); % this is right
vvad = -speed .* cosd(direction);
J = load_edson_data_10min();
umast = -J.wspd .* sind(J.wdir); % this is also right, and they agree
vmast = -J.wspd .* cosd(J.wdir);
uship =  J.SOG  .* sind(J.COG ); % course in direction of motion, so no - sign
vship =  J.SOG  .* cosd(J.COG );
height_mid_wind = 0.5*( height_wind(1:end-1) + height_wind(2:end) );
dz_wind = diff(height_wind);
dz_smth = dz_wind;
dz_smth(2:end) =  (height_mid_wind(2:end)/2.2e-3).^(1/4.2);
shear2 = (diff(uvad,1)./dz_smth).^2 + (diff(vvad,1)./dz_smth).^2;

ii = isfinite(J.yday);
% relative wind profiles
urel = uvad' - interp1(J.yday(ii), uship(ii), yDay_wind(:), 'linear');
vrel = vvad' - interp1(J.yday(ii), vship(ii), yDay_wind(:), 'linear');
% ship maneuvers not removed!!
% horizontal ship-relative speed [m/s] from lidar VADs and interpolated ship velocity
urel_wind = sqrt(urel.^2 + vrel.^2);

%% allocate output data variables
nhr = 24 * 6;
n10min = nhr * 3;
kmax = 60; % find(zData<2e3,1,'last')
kradarmax = 66;
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

%% loop over hourly files
% handle vertical stares split across hourly files
itmax   = 4000;
itrmax  = 12700;
yday    = NaN(itmax, 1);
W       = NaN(itmax, kmax);
ydradar = NaN(itrmax, 1);
Z       = NaN(itrmax, kradarmax);

% initialize trailing previously-loaded part of chunk
lastyday      = NaN(600, 1);
lastW         = NaN(600, kmax);
lastydayradar = NaN(1000, 1);
lastZ         = NaN(1000, kradarmax);
[ydstart, ydend] = deal(zeros(432,1));
ntlast      = 0; % number of time points to concatenate from prev hour
ntlastradar = 0;
bigind      = 0; % counter for filling output arrays
% must zero ntlast, bigind before starting ifile loop!

for ifile = istart : length(files)
    fprintf(1,'%s\n',files(ifile).name)
    retrfile  = fullfile(way_DopplerData,files(ifile).name);
    yDay      = ncread(retrfile,'yDay'); % not same as horizontal wind height
    mcVelData = ncread(retrfile,'mcVelData',[1 1],[kmax Inf])';
    snrData   = ncread(retrfile,'snrData',[1 1],[kmax Inf])';
    mcVelData(snrData<-7) = NaN; % filter out bogus noisy velocities
    
    % read radar to exclude (include?) clouds from dissipation calculation
    ix = min(100, length(yDay));
    yyday = floor( yDay(ix) );
    hhour = floor(mod((yDay(ix) - yyday)*24, 24));
    radarfilestx = dir(fullfile(way_RadarData,sprintf('2011%03i%02i*MMCRMom.nc',yyday,hhour)));
    
    if ~isempty( radarfilestx )
        radcount = 0;
        ydayradar = NaN(radcount); Zradar = NaN(radcount, kradarmax); % minimalist preallocation
        for iradarfile = 1 : length(radarfilestx) % loop through radar files in the hour
            radarfile   = radarfilestx(iradarfile).name;
            base_time   = ncread(fullfile(way_RadarData, radarfile), 'base_time'  );
            time_offset = ncread(fullfile(way_RadarData, radarfile), 'time_offset');
            tt  = time_offset + double(base_time);
            nta = length(tt); % time length of current file
            ydayradar(radcount+(1:nta),1) = tt / (60 * 60 * 24) + datenum(1970, 1, 1) - datenum(year, 0, 0);
            Zradar   (radcount+(1:nta),:) = double( ncread( [way_RadarData '/' radarfile], 'Reflectivity', [1 1], [kradarmax Inf] ))';
            radcount = radcount + nta;
        end % loop through radar files in the hour
        
        % reconstruct strictly increasing time ydayradar
        iendcont   = [find(abs(diff(ydayradar))>1); length(ydayradar)]; % end of contiguous regions
        istartcont = [1; iendcont(1:end-1)+1];
        for i = 1 : length(iendcont)
            ii  = istartcont(i):iendcont(i);
            nii = length(ii);
            ydayradar(ii) = linspace(ydayradar(istartcont(i)), ydayradar(iendcont(i)), nii);
        end
        hourradar = ( ydayradar - fix(ydayradar(1)) ) * 24; % subtract off JD and convert to an hour
        
        % concatenate data onto previous trailing chunk
        ntnext = length(yDay);
        nthour = ntlast + ntnext;
        yday(1:nthour)   = [lastyday(1:ntlast); yDay];
        W   (1:nthour,:) = [lastW(1:ntlast,:); mcVelData(:,1:kmax)];
        
        ntnextradar = length(ydayradar);
        nthourradar = ntlastradar + ntnextradar;
        ydradar(1:nthourradar) = [lastydayradar(1:ntlastradar); ydayradar];
        Z(1:nthourradar,:) = [lastZ(1:ntlastradar,:); Zradar];
        
        % exclude clouds in the HRDL Doppler vel data
        ii = [true; diff( ydradar(1:nthourradar) ) > 0]; % interpolate across when radar time goes backwards
        iscloud = bsxfun(@gt, Z(ii,:), noisefloor')';
        [X,Y]   = meshgrid(ydradar(ii), hradar);
        [X1,Y1] = meshgrid(yday(1:nthour),zData(1:kmax));
        iscloudlidar = (interp2(X,Y,double(iscloud), X1,Y1, 'nearest')>0)';
        [i,j] = ind2sub(size(iscloudlidar), find(iscloudlidar));
        isabovecloud = iscloudlidar;
        uj = unique(j);
        for pj = 1 : length(uj)
            isabovecloud( min(i( j==uj(pj) )) : kmax,uj(pj) ) = true;
        end
        % exclude HRDL data above any clouds, since lidar will be attenuated
        W(isabovecloud) = NaN;
        
        % loop through, interp and calc spectra for 10-min chunks
        %starts=[0; find(diff(yday)>0.005)]+10; % find, chop off start of each 10min chunk
        %ends=[find(diff(yday)>0.005)-1; length(yday)];
        %ends=ends(1:length(starts));
        % find chunk start and end indices from diff(yday) gaps of more
        % than 5 minutes
        dift    = find( diff(yday(1:nthour)) > 5/(24*60) );
        starts  = [1; dift+10]; % chop off 1st 10 points of each vertical stare
        starts(diff(starts) <= 10) = []; % remove corner case of contiguous diff(yday) gaps
        ends    = dift - 1;
        nchunks = length(ends);
        startsradar = floor(interp1( ydradar(1:nthourradar), 1:nthourradar, yday(starts) ));
        endsradar   = ceil (interp1( ydradar(1:nthourradar), 1:nthourradar, yday(ends)   ));
        
        % loop through chunks and do analysis
        for ind = 1 : nchunks % index of the 10 minute interval in the hour
            i10 = starts(ind) : ends(ind); % first index in the 10 minute interval
            nt  = length(i10);
            
            if nt >= di / 4 % skip to next iteration if window too small
                bigind          = bigind + 1;          % increment index
                ydstart(bigind) = yday( starts(ind) );
                ydend(bigind)   = yday(   ends(ind) );
                ydayw10(bigind) = yday( starts(ind) ); % current time of this 10-minute window
                nw              = sum( isfinite(W(i10,:)) );
                [nmax,knmax]    = max( nw ); % nw is the max number of nonnoise returns, knmax is its height index
                
                % choose the highest vertical level ik that has >0.25*nmax of this coverage
                ik  = min(kmax, find(nw>0.25*nmax, 1, 'last'));  % highest
                hk  = find(nw>0.5*nmax, 1, 'first'); % lowest
                kk  = (max(hk-1, 1) : min(ik+1, kmax))';
                nk  = length(kk);
                w10 = W(i10, kk);
                
                % filter out extreme values
                w10( bsxfun(@gt, abs(bsxfun(@plus,w10,-nanmean(w10))), 4*nanstd(w10)) ) = NaN;
                iwt = isfinite(w10);
                count(bigind,kk) = sum(iwt);
                
                % interpolate horizontal wind to the middle of each chunk
                ydmid = nanmedian(yday(i10));
                iiw   = isfinite(yDay_wind);
                jjw   = isfinite(height_wind);
                Ulev  = interp2(yDay_wind(iiw), height_wind(jjw), urel_wind(iiw,jjw)', ydmid,zData(1:kmax), 'linear');
                
                for k = 1 : nk % loop through vertical levels
                    if count(bigind,kk(k)) > di / 4
                        S(:) = fast_acov_spectr(w10(:,k), nwin); % Ignore Nyquist-frequency power estimate!
                        if isfinite(S(1))
                            spec(bigind,:,kk(k)) = S;
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
                            sigv2(bigind,kk(k)) = dF * sum( S(1:knoise(bigind,kk(k))) - Snoise(bigind,kk(k)) );
                            ulev (bigind,kk(k)) = Ulev(k);
                        end % finite spectra were returned
                    end % there is data for spectra
                end % k level where there is data
                
            end % there is 10 min HRDLvelocity data
        end % 10 min
        
        % save last trailing chunk fragment to concatenate to next hour
        ilast  = starts(end) : nthour;
        ntlast = length(ilast);
        lastyday(1:ntlast)   = yday(ilast);
        lastW   (1:ntlast,:) = W   (ilast,:);
        % trailing chunk indices for radar
        if startsradar(end) > 0 % radar file goes to end of hour
            ilastradar  = startsradar(end) : nthourradar;
            ntlastradar = length(ilastradar);
            lastydayradar(1:ntlastradar)   = ydradar(ilastradar);
            lastZ        (1:ntlastradar,:) = Z      (ilastradar,:);
        else
            ntlastradar = 0;
            lastydayradar = NaN;
            lastZ(1,:)    = NaN;
        end
    end % at least one radar file exists block
end     % HRDL  file

% truncate output data
ydayw10 = ydayw10(1:bigind);
ydstart = ydstart(1:bigind);
ydend   = ydend  (1:bigind);
epsilon = real(epsilon);
epsilon = epsilon(1:bigind,:);
sigv2   = sigv2  (1:bigind,:);
count   = count  (1:bigind,:);
spec    = spec   (1:bigind,:,:);
knoise  = knoise (1:bigind,:);
Snoise  = Snoise (1:bigind,:);
Sthresh = Sthresh(1:bigind,:);
ulev    = ulev   (1:bigind,:);

%% save data
% save hrdlturb4_512n_alldays.mat ydayw10 ydstart ydend F spec epsilon sigv2 count ulev knoise Snoise Sthresh zData

%% plot nondimensional spectra scaled by Kolmogorov power spectral density (dissipation*dynvisc^5)^1/4
% for cloud top
dynvisc = 1.63e-5; %m^2/s, at ~900 hPa
Kpsd = sqrt(sqrt( epsilon *  (dynvisc^5) ));  % Kolmogorov power spectral density (dissipation*dynvisc^5)^1/4  [U^2 / L^-1]
Kks  = sqrt(sqrt( epsilon ./ (dynvisc^3) ));  % Kolmogorov wavenumber (dissipation*dynvisc^3)^1/4  [L^-1]
eta = 1 ./ Kks;                            % Kolmogorov length scale
ulevx = reshape(ulev,[size(ulev,1) 1 size(ulev,2)]);
Fx = reshape(repmat(F,[1 size(ulevx,1)])',[size(ulevx,1),length(F),1]);
wavenumber = 2 * pi * Fx ./ ulevx;
kscaled  = wavenumber ./ reshape(Kks,[size(ulev,1) 1 size(ulev,2)]); % nondimensional wavenumber is scaled by Kolmogorov wavenumber
speck    = spec     .* ulevx/(2*pi);   % wavenumber spectrum (spec is frequency spectrum!)
Sfscaled = spec     ./ reshape(Kpsd,[size(ulev,1) 1 size(ulev,2)]) ;  % nondimensional (TKE) power frequency-spectral density scaled by Kolm'v PSD
Skscaled = Sfscaled .* ulevx/(2*pi);  % nondimensional (TKE) power wavenumber-spectral density scaled by Kolm'v PSD

% cut off noise
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

% plot time spectrogram
lev=5;
clf
pcolor(log10(Skscaled(:,:,lev))'); shading flat; set(gca,'yscale','log')
%pcolor(log10(Scut(:,:,lev))'); shading flat; set(gca,'yscale','log')
hold on
plot(knoise(:,lev),'r.-')

% plot vertical spectrogram
tim = 105; % 2011 Oct 09
Fnoise = zeros(60,0);
Fnoise(knoise(tim,:)>0) = F(knoise(tim,knoise(tim,:)>0));
clf
pcolor(1:60, F, squeeze( log10(Skscaled(tim,:,:)) ) ); shading flat; set(gca,'yscale','log')
% pcolor(squeeze( log10(Scut(tim,:,:)) ) ); shading flat; set(gca,'yscale','log')
ylabel('frequency index'); xlabel('veritcal level')
hold on
plot(Fnoise,'r.-')
plot(1e3*epsilon(tim,:),'bo')

tim = 105;
clf
loglog(squeeze(Skscaled(tim,:,2:2:8)),'b')
hold on
plot(knoise(tim,2:2:8),squeeze(Skscaled(tim,knoise(tim,2:2:8),2:2:8)),'bo')
loglog(squeeze(Skscaled(tim,:,45:46)),'r')
plot(knoise(tim,45:46),squeeze(Skscaled(tim,knoise(tim,45:46),45:46)),'ro')

% top returns are flat at low wavenumber: don't look like turbulence

% plot scaled not by energy containing scale or by noise, but by inferred
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
title('HRDL DYNAMO 750 m Kolmogorov-scaled w wavenumber spectra', 'fontweight','normal')
ylabel('power spectral density   S_{ww}(k)/(\epsilon\nu^5)^{1/4}')
xlabel('wavenumber   k\eta')

% plot compensated spectrum
% clf
loglog(kscaled(:,:,lev)',Skscaled(:,:,lev)'.*kscaled(:,:,lev).^(5/3)','.','markersize',1,'color',0.7*[1 1 1])
hold on
loglog(kscaled(:,2:end,lev)',Scut(:,2:end,lev)'.*kscaled(:,2:end,lev).^(5/3)','b.','markersize',1)
plot([1e-6 1e-2],C1prime*[1 1],'r-')
% saveas(gcf,'Kolmogorov_spectrum4all.eps','epsc')

%{
errorbar(nanmean(Scut(:,4:end)'.*kscaled(:,4:end).^(5/3)',2),nanstd(Scut(:,4:end)'.*kscaled(:,4:end)'.^(5/3),2)./sqrt(sum(isfinite(Scut(:,4:end)'.*kscaled(:,4:end)'),2)))
% k,kscaled changes with the wind, so binavg
kbin=logspace(5e-5,5e-2,31);
[mb,sb,wb,nb]=binavg(kbin,kscaled(:,4:end),Scut(:,4:end).*kscaled(:,4:end).^(5/3));
[mlb,slb,wlb,nlb]=binavg(kbin,kscaled(:,4:end),log(Scut(:,4:end).*kscaled(:,4:end).^(5/3)));
errorbar(kbin,mb,sb./sqrt(nb))
set(gca,'xscale','log')
%}

% diagnose vertical noise in wind profiles
%{
hold on
plot(dz_wind(kwind), height_mid_wind(kwind), '.-')

clf
loglog(dz_wind(kwind), height_mid_wind(kwind), '.-')
hold on
loglog([9 30], 2.2e-3*[9 30].^(4.2))
loglog(height_mid_wind(kwind), '.-')
%}

% plot dissipation of the mixed layer
i1=ydayw10<305; i2=ydayw10>305;
b2rcolormap(16);
clf
ax(1) = subplot(2,1,1, 'align');
pcolor(ydayw10(i1)-0.0069,zData(1:kmax)/1e3,log10(epsilon(i1,:)')); shading flat
hold on
%contour(ydayw10(i1),zData(1:kmax)/1e3,count(i1,:)',[200 400 800],'k'); shading flat
plot(J.yday, -J.Solardn/1e3,'k-','linewidth',1.5)
plot(J.yday,  J.wspd/10,'b-','linewidth',1)
plot(J.yday, (J.SST-28)/2,'r-','linewidth',1)
caxis([-6 -3]); ylim([0 2])
colorbar
ylabel('height (km)')
set(gca, 'color',0.7*[1 1 1], 'fontsize',14, 'tickdir','out')
ax(2) = subplot(2,1,2, 'align');
pcolor(ydayw10(i2)-0.0069,zData(1:kmax)/1e3,log10(epsilon(i2,:)')); shading flat
hold on
%contour(ydayw10(i2),zData(1:kmax)/1e3,count(i2,:)',[200 400 800],'k'); shading flat
plot(J.yday, -J.Solardn/1e3,'k-','linewidth',1.5)
plot(J.yday,  J.wspd/10,'b-','linewidth',1)
plot(J.yday, (J.SST-28)/2,'r-','linewidth',1) % sst ylim [28 32]
caxis([-6 -3]); ylim([0 2])
colorbar
xlabel('yearday 2011')
ylabel({'solar flux (kW m^{-2})'; 'wind speed (dm s^{-1})'})
set(gca, 'color',0.7*[1 1 1], 'fontsize',14, 'tickdir','out')
orient landscape
% saveas(gcf,'dissipation_ML_512n.eps','epsc')
% saveas(gcf,'dissipation_ML_512n.png','png')

% show diurnal warm layers when wind is weak
% plot shear^2 in the mixed layer
iwind = yDay_wind >= 316 & yDay_wind <= 322;
kwind = height_mid_wind <= 2200;
clf
ax(1) = subplot(2,1,1, 'align');
contourf(yDay_wind(iwind), height_mid_wind(kwind)/1e3, log10(shear2(kwind, iwind)), -9:0.5:-2,'edgecolor','none')
axis([316 322 0 2])
caxis([-7 -2])
ylabel({'shear^2 (m^2 s^{-1})', 'height (km)'}); %xlabel('2011 yearday');
set(gca,'fontsize',14, 'tickdir','out', 'xaxislocation','top')
ax(1).Position(2) = ax(1).Position(2) - 0.08;
colorbar
text(316,2.5, 'Nov 13', 'clipping','off', 'horizontalalignment','center', 'fontsize',14)
text(319,2.5, 'Nov 16', 'clipping','off', 'horizontalalignment','center', 'fontsize',14)
text(322,2.5, 'Nov 19', 'clipping','off', 'horizontalalignment','center', 'fontsize',14)
% saveas(gcf, 's2.png')
ax(2) = subplot(2,1,2, 'align');
pcolor(ydayw10(i2)-0.0069,zData(1:kmax)/1e3,log10(epsilon(i2,:)')); shading flat
hold on
% contour(yDay_wind(iwind), height_mid_wind(kwind)/1e3, log10(shear2(kwind, iwind)), [-3 -3],'edgecolor','y')
%contour(ydayw10(i2),zData(1:kmax)/1e3,count(i2,:)',[200 400 800],'k'); shading flat
plot(J.yday, -J.Solardn/1e3,'k-','linewidth',1.5)
plot(J.yday,  J.wspd/10,'b-','linewidth',1)
plot(J.yday, (J.SST-28)/2,'r-','linewidth',1) % SST  ylim [28 32]
plot(J.yday, (J.T10-27)/2,'g-','linewidth',1) % Tair ylim [27 31]
caxis([-6 -3]); ylim([0 2])
cb(2) = colorbar; cb(2).YLabel.String = 'log_{10} dissipation';
xlim([316 322]) % Nov 13-19
xlabel('yearday 2011')
ylabel({'solar flux (kW m^{-2})'; 'wind speed (dm s^{-1})'})
set(gca, 'color',0.7*[1 1 1], 'fontsize',14, 'tickdir','out')
orient landscape

% other diagnostics follow
%{
% retrievals with fewer than 800 counts are suspect.
ii = count<800 & epsilon>10^-4 & zData(1:kmax,ones(length(count),1))'>1e3;
[it,iz] = ind2sub(size(count),find(ii));
x = zeros(nwin/2, length(it));
for i = 1:length(it)
    x(:,i) = spec(it(i),:,iz(i));
end
% loglog(x);

for i = 1:length(it); y(i) = count(it(i),iz(i)); end

mask = NaN(size(epsilon));
for i = 1:length(it)
    mask(it(i),iz(i)) = -6;
end
h = pcolor(ax(1),ydayw10(i1)-0.0069,zData(1:kmax)/1e3,mask(i1,:)'); h.EdgeColor='none';
h = pcolor(ax(2),ydayw10(i2)-0.0069,zData(1:kmax)/1e3,mask(i2,:)'); h.EdgeColor='none';
% covers lots of bad cases near top of retrievals

% most suspect
farr = find(x(1,:)>1e3); % 25
% I think shows the spectral response of the windowing, redder than other
% spectra.
y = zeros(nwin/2, length(farr));
for i = 1:length(farr)
    y(:,i) = spec(it(farr(i)),:,iz(farr(i)));
end
% not referencing back to epsilon(time, height) properly !?!?
plot(ax(1),ydayw10(it(farr)), zData(iz(farr))/1e3, 'gx','markersize',8);
plot(ax(2),ydayw10(it(farr)), zData(iz(farr))/1e3, 'gx','markersize',8);

% plot spectrogram in order of number of good data
[sortcount, ordercount] = sort(count(:), 'ascend');
isf = isfinite(epsilon(ordercount));
ocf = ordercount(isf);
[ot, oh] = ind2sub(size(count),ocf);
[speco, speko, ksco] = deal( zeros(nwin/2, length(ot)) );
for i=1:length(ot) % loops are comprehensions, so easy!
    speco(:,i) =     spec(ot(i),:,oh(i))';
    speko(:,i) = Skscaled(ot(i),:,oh(i))';  % Kolmogorov scaled wavenumber spectra ordered by count
    ksco (:,i) =  kscaled(ot(i),:,oh(i))';  % Kolmogorov scaled wavenumber
end

clf; ax=[];
ax(1) = subplot(2,1,1);
pcolor(sortcount(isf), F     , log10(speco)); shading flat; set(gca,'yscale','log','fontsize',14)
ylabel('frequency [s^{-1}]'); xlabel('count in window')
ax(2) = subplot(2,1,2);
% pcolor(1:sum(isf), ksco, log10(speko));
% mesh(1:sum(isf), ksco, log10(speko));
% waterfall(1:sum(isf), ksco, log10(speko)); % crashes
shading flat; set(gca,'yscale','log','fontsize',14)
% set(gca, 'xtick',1000:1000:10000, 'xticklabel',count(ocf(1000:1000:10000)))
set(gca,'color',[0.5 0.5 0.5])
caxis([5 8]); ylim([2e-5 2e-2])
ylabel('nondimensional wavenumber'); xlabel('count in window')
hold on
plot(sortcount(isf)/1e5)
orient landscape
% saveas(gcf, 'spectrogram_vs_counts.png')

% plot counts
clf
subplot(2,1,1)
pcolor(ydayw10(i1),zData(1:kmax)/1e3,count(i1,:)'); shading flat
colorbar
ylabel('height (km)')
set(gca,'color',0.7*[1 1 1])
subplot(2,1,2)
pcolor(ydayw10(i2),zData(1:kmax)/1e3,count(i2,:)'); shading flat
colorbar
xlabel('yearday 2011')
set(gca,'color',0.7*[1 1 1])


% what's wrong with the dissipation at the highest few bins?
axis([218.5 282.3 1.3 1.9])

[i,j] = ind2sub(size(count),find(count<600 & count>25));
x = permute( kscaled,[2 1 3]);
y = permute(Skscaled,[2 1 3]);
ii = count<700 & count>25;
loglog(x(:,ii),y(:,ii),'.')
hold on
loglog(x(:,~ii),y(:,~ii),'c.')
%}


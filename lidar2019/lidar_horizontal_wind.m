% Compute horizontal wind from Streamline XR User5_06 scans.
% Doppler velocity convention is positive towards the lidar.
switch computer
    case 'PCWIN64'
        addpath(genpath('C:\Users\Simon de Szoeke\Documents\MATLAB\user_tools\graphics\'));
    case 'MACI64'
        addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/graphics/'));
end

nrange = 104; % number of range gates
rangegate = 30.0;
kts2ms = 0.51444444; % factor converts knots to m/s

pltit = false;
prtit = true;

nbig = 6*24*20;
[uwnd, vwnd, uerr, verr, urel, vrel] = deal( NaN(nrange, nbig) );
shear2 = deal( NaN(nrange-1, nbig) );
[time, muship, mvship, mheading, mcog, msog] = deal( NaN(nbig,1) );

% files to read
% metdir = '/Volumes/cruise/SR1911/share/data/MET/';
metdir = '~/Data/cruises/MISOBOB_2019/SR1911/ship/met/mat/';
datadir = '~/Data/cruises/MISOBOB_2019/SR1911/lidar/HaloLidar_sn0610-06_leg2/';
dd = dir(fullfile( datadir, '2019*' ));
dates = {dd.name};

ic = 1; % start index

% file = '/Volumes/MISOBOB/Halo Lidar sn0610-06 LEG 2/20190712/Stare_06_20190712_00.hpl';
for id = 1:length(dates) % start on 20190708
    fprintf(1, '%s', dates{id})
    met = getfield( load( fullfile(metdir, ['MET_' dates{id}(3:end) '.mat']) ), 'MET');
    
    for hour = 0:23
        hh = sprintf('%02i', hour);
        fprintf(1, '|')
        % starefile = fullfile(datadir, dates{id}, ['Stare_06_' dates{id} '_' hh '.hpl']);
        scanfiles  = dir( fullfile(datadir, dates{id}, ['User5_06_' dates{id} '_' hh '*.hpl']) );
        mday0 = datenum( dates{id}, 'yyyymmdd' );
        
        for isc = 1:length(scanfiles)
            fprintf(1, '.')
            % read the stare file into the structure a
            file = fullfile(scanfiles(isc).folder, scanfiles(isc).name);
            mm = scanfiles(isc).name(end-7:end-6);
            b = read_streamlinexr(file, nrange);
            
            % mask out velocities where the SNR is too low
            thr = 1.0; % lowered threshold to get more mean wind data
            mask = b.intensity > thr; % SNR threshold
            nx = sum(mask,2);
            nanmask = zeros(size(mask));
            nanmask(~mask) = NaN;
            
            mday = mday0 + b.dechour / 24;
            time(ic) = nanmean(mday);
            
            % get ship heading, speed, and course over ground at ray times
            [tu iu] = unique( met.time );
            heading =    interp1( tu, met.GY(iu), mday );
            cog =        interp1( tu, met.CR(iu), mday );
            sog = kts2ms*interp1( tu, met.SP(iu), mday );
            uship   = sog .* sind( cog );  % to add to urel,vrel
            vship   = sog .* cosd( cog );
            
            % mean quantities to be saved
            mheading(ic) = nanmean(heading);
            mcog(ic)     = nanmean(mcog);
            msog(ic)     = nanmean(sog);
            muship(ic)   = nanmean(uship);
            mvship(ic)   = nanmean(vship);
            
            % trig of azimuth and elevation angles
            cosel = cosd( b.elevation  );
            phi = b.azimuth - 90.0 + heading;
            % b.azimuth = 0 points to port; +90 degrees to the bow.
            % Adjusted by -90 degrees so phi = 0 points North.
            x = cosd( phi ) .* cosel; % northing
            y = sind( phi ) .* cosel; % easting
            
            % find urel, vrel relative to lidar in earth coordinate systen
            % negative dopplervel switches from neg-away convention to positive-away
            % urel > 0 to starboard
            % vrel > 0 forward
            xg = bsxfun(@times, mask, x');
            sxx = nansum( xg .* xg, 2 );
            yg = bsxfun(@times, mask, y');
            syy = nansum( yg .* yg, 2 );
            dvel = mask.*(b.dopplervel); % (-b.dopplervel); ???
            % assumes zero divergence:
            vrel(:,ic) = ( dvel * x ) ./ sxx; % northing
            urel(:,ic) = ( dvel * y ) ./ syy; % easting
            %       inner prod. sums
            
            % maybe should allow a nonzero intercept and compute divergence
            % The mean divergence and deformation might be assumed to be 
            % constant with height, so can be solved from du/dx, dv/dy.
            % Shear and vorticity cannot be estimated from radial
            % velocities.
            % 
            % rxy = cosd(b.elevation) * 30.0*(0:(nrange-1)): % horizontal displacement, m
            
            % error estimates of slope for goodness of fit
            verr(:,ic) = sqrt( nansum(mask.*(dvel - vrel(:,ic)*y').^2, 2) ./ ( (nx-2) .* (sxx - sum(xg,2).^2./nx) ) );
            uerr(:,ic) = sqrt( nansum(mask.*(dvel - urel(:,ic)*x').^2, 2) ./ ( (nx-2) .* (syy - sum(yg,2).^2./nx) ) );
            
            % add ship velocity
            uwnd(:,ic) = urel(:,ic) + nanmean(uship); % not done: ship vel should be added to each ray individually if it is manuevering
            vwnd(:,ic) = vrel(:,ic) + nanmean(vship);
            
            du = diff(urel(:,ic));
            dv = diff(vrel(:,ic));
            shear2(:,ic) = (du.*du + dv.*dv) / (rangegate * mean(sind( b.elevation ))).^2; % s^-1
            
            %% plot the 10-min wind
            height_wnd = 1e-3*rangegate*(0:nrange-1) * mean(sind( b.elevation )); % km
            
            if pltit % very slow
                clf; clear ax;
                ax(1) = subplot(2,2,1);
                plot(uwnd(:,ic),vwnd(:,ic),'.-'); xlabel('u (m/s)'); ylabel('v (m/s)');
                hold on
                plot(uwnd(1:5:end,ic), vwnd(1:5:end,ic), 'o');
                title('lidar wind profile')
                ax(2) = subplot(2,2,3);
                plot([0 0], [0 3],'k-')
                hold on
                plot(uwnd(:,ic), height_wnd, '.-')
                plot(uwnd(1:5:end,ic), height_wnd(1:5:end), 'o')
                plot(5e3*shear2(:,ic), (height_wnd(1:end-1)+height_wnd(2:end))/2, '-','linewidth',2)
                xlim([-20 20])
                xlabel('u (m/s)')
                ylabel('height (km)')
                ax(3) = subplot(2,2,4);
                plot([0 0], [0 3],'k-')
                hold on
                plot(vwnd(:,ic), height_wnd, '.-')
                plot(vwnd(1:5:end,ic), height_wnd(1:5:end), 'o')
                xlim([-20 20])
                xlabel('v (m/s)')
                set(ax(:), 'fontsize',12)
                if prtit
                    print('-depsc', ['./plots/lidar_wind_' dates{id} '-' hh mm '.eps']);
                end
            end
            
            ic = ic + 1; % increment index
            
        end
        
    end
    fprintf(1,'\n')
end

% truncate
ii = 1:ic-1;
time = time(ii);
muship = muship(ii);
mvship = mvship(ii);
mheading = mheading(ii);
mcog = mcog(ii);
msog = msog(ii);
urel = urel(:,ii);
vrel = vrel(:,ii);
uwnd = uwnd(:,ii);
vwnd = vwnd(:,ii);
uerr = uerr(:,ii);
verr = verr(:,ii);

save('lidar_horizontal_wind.mat', 'time', 'height_wnd', 'urel', 'vrel', 'uwnd', 'vwnd', 'uerr', 'verr', 'muship', 'mvship', 'mheading', 'mcog', 'msog','shear2')

%% this is how I would interpolate to the stares:
% wnd = load('/Users/sdeszoek/Data/cruises/MISOBOB_2019/SR1911/scientific_analysis/programs/lidar/lidar_horizontal_wind.mat');
% interp2( wnd.time, wnd.height_wnd, wnd.urel, wnd.vrel, stare.time );

%% plot
b2rcolormap(21);
clf
subplot(2,1,1)
pcolor(time-datenum(2019,7,0), height_wnd(1:34), uwnd(1:34,:)); shading flat;
caxis([-20 20])
colorbar
hold on 
plot(met.time-datenum(2019,7,0), -met.TW.*sind(met.TI)/20,'k-')
plot(time-datenum(2019,7,0), uwnd(5,:)/20,'y')
colorbar
subplot(2,1,2)
pcolor(time-datenum(2019,7,0), height_wnd(1:34), vwnd(1:34,:)); shading flat;
caxis([-20 20])
colorbar
hold on 
plot(met.time-datenum(2019,7,0), -met.TW.*cosd(met.TI)/20,'k-')
plot(time-datenum(2019,7,0), vwnd(5,:)/20,'r')


pcolor(time-datenum(2019,7,0), height_wnd(1:34), vrel(1:34,:)); shading flat;


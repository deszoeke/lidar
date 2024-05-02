% run lidar_horizontal_wind first --> lidar_horizontal_wind.mat

% cd('/Users/sdeszoek/Data/cruises/MISOBOB_2019/SR1911/scientific_analysis/programs/lidar'0
paren = @(X, varargin) X(varargin{:});

switch computer
    case 'PCWIN64'
        cruiseshare = '\\sr-nas1.rv-sallyride.ucsd.edu\cruise\';
        metdir = fullfile( cruiseshare, 'SR1911/share/data/MET/' );
        datadir = 'D:/Halo Lidar sn0610-06 LEG 2/';

        addpath(genpath('C:\Users\Simon de Szoeke\Documents\MATLAB\user_tools\graphics\'));
        cd('D:\programs\lidar')
    case 'MACI64'
        % cruiseshare = '/Volumes/cruise/'; % not available
        % datadir = '/Volumes/MISOBOB/Halo Lidar sn0610-06 LEG 2/';
        metdir = '/Users/sdeszoek/Data/cruises/MISOBOB_2019/SR1911/ship/met/mat/';
        datadir = '/Users/sdeszoek/Data/cruises/MISOBOB_2019/SR1911/lidar/HaloLidar_sn0610-06_leg2/';
        addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/graphics/'));
        cd('/Users/sdeszoek/Data/cruises/MISOBOB_2019/SR1911/scientific_analysis/programs/lidar');
end
addpath('./turb_subs/')

% metdir = fullfile(cruiseshare, '/SR1911/share/data/MET/');

plotit = false & true;

nrange = 104; % number of range gates
rangegate = 30.0;
range = 1e-3*rangegate*(0:(nrange-1)'); % km from lens
height = range + 0.015; % height of vertical stares, km
fs = 2.0;  % sampling frequency 1/s
thr = 1.5; % noise threshold seems to work from experimentation

kts2ms = 0.51444444; % factor converts knots to m/s

dd = dir(fullfile( datadir, '2019*' )); % data folders by day
dates = {dd.name};

% load cruise lidar horizontal wind data
% H = load('D:/programs/lidar/lidar_horizontal_wind.mat');
H = load('/Users/sdeszoek/Data/cruises/MISOBOB_2019/SR1911/scientific_analysis/programs/lidar/lidar_horizontal_wind.mat');

% load sounding horizontal wind data
% sondedir = '/Volumes/cruise/SR1911/share/data/radiosonde/gridded/';
% sondedir = '\\sr-nas1.rv-sallyride.ucsd.edu\cruise\SR1911\share\data\radiosonde\gridded\';
sondedir = '/Users/sdeszoek/Data/cruises/MISOBOB_2019/SR1911/radiosonde/gridded/';
P = getfield( load(fullfile(sondedir, 'SallyRide_grd.mat')), 'grd');

nbig = (length(dates)-2)*24*6;
[epsilon, Ulev] = deal( NaN( nrange, nbig ) ); % allocate epsilon
[ time_start, time_end ] = deal( NaN( nbig, 1 ) );

count = 0;
% example
% file = '/Volumes/MISOBOB/Halo Lidar sn0610-06 LEG 2/20190712/Stare_06_20190712_00.hpl';
for id = 11:length(dates) % daily folders, start on 20190716
    fprintf(1, '%s', dates{id})
    flist = dir(fullfile( dd(id).folder, dd(id).name, 'Stare_06*.hpl' ));
    % get met file for ship heading
    met = getfield( load( fullfile(metdir, ['MET_' dates{id}(3:end) '.mat']) ), 'MET');

    for ifile = 1:length(flist) % hourly files
        fprintf(1, '.%s', flist(ifile).name(19:20))
        file = fullfile( flist(ifile).folder, flist(ifile).name );
        S = read_streamlinexr(file, nrange); % slow
        mtime = datenum(flist(ifile).name(end-14:end-4), 'yyyymmdd_HH') + mod(S.dechour,1)/24;
        % just one hour
        
        % roll and pitch are reversed
        % Streamline XR is rotated on deck.
        %    roll vector is to port, (right-handed) bow  tilts down
        %   pitch vector is to bow , (right-handed) port tilts up 
        
        % interp relative horiz wind in earth coordinates to lidar stare grid
        urel = interp2( H.time, H.height_wnd, H.urel, mtime, range );
        vrel = interp2( H.time, H.height_wnd, H.vrel, mtime, range );
        
        % interp sounding true wind to lidar stare (vertical displacement) grid
        nhmax = 150;
        ntsnd = size(P.h, 2);
        hint = nanmean(P.h,2); % height interp. vector
        ii = isfinite(hint);
        ui = interp1( hint(ii), P.u(ii,:), 1e3*height, 'linear' );
        vi = interp1( hint(ii), P.v(ii,:), 1e3*height, 'linear' );
        usnd = interp1( P.time, ui', nanmean(mtime), 'linear' )';
        vsnd = interp1( P.time, vi', nanmean(mtime), 'linear' )';
        
        % ship velocity
        % get ship heading, speed, and course over ground at stare times
        [tu, iu] = unique( met.time );
        heading =    interp1( tu, met.GY(iu), mtime );
        cog =        interp1( tu, met.CR(iu), mtime );
        sog = kts2ms*interp1( tu, met.SP(iu), mtime );
        uship   = sog .* sind( cog );  % to add to urel,vrel
        vship   = sog .* cosd( cog );

        % ship relative velocities from soundings
        urelsnd = usnd - uship';
        vrelsnd = vsnd - vship';
        Urelsnd = sqrt( urelsnd.*urelsnd + vrelsnd.*vrelsnd );
        
        % rotate relative velocity from earth- to ship-coordinate system
        [ut, ui] = unique(met.time);
        hd = interp1(ut, met.GY(ui), mtime); % heading
        sh = sind(hd); ch = cosd(hd);
        urelshp =  urel.*ch' + vrel.*sh'; % to starboard
        vrelshp = -urel.*sh' + vrel.*ch'; % to bow   % range, time
        Urelshp = sqrt( urelshp.*urelshp + vrelshp.*vrelshp );
        
        % resolve relative wind into Doppler vel using Streamline XR roll and pitch
        rayvel = vrelshp.*sind(S.roll)' + urelshp.*sind(S.pitch)';
        % mean wind in radial (ray) direction, positive away from lidar

        % mask the Doppler velocities
        mask = S.intensity > thr;
        nanmask = zeros( size(mask) );
        nanmask(~mask) = NaN;
        Dvm = nanmask + S.dopplervel; % masked velocity array

        subtractmeanwindcomponent = false;
        if subtractmeanwindcomponent % subtract mean relative wind from Doppler velocity
            x = Dvm - rayvel;
        else
            x = Dvm;
        end
        % subtract ray-mean velocity
        Dv = x - nanmean(x, 1);
        % Have not subtracted lidar motion (mostly heave) explicitly.
        % VectorNav stable table data is in
        % /Volumes/MISOBOB/Halo Lidar sn0610-06 LEG 2/Stable_Data/
        % I think it suffices to subtract the ray-mean velocity.
        
        % divide into continuous 7.5 min chunks of 2 Hz data
        chlast  = [find(diff(mtime) > 1e-3); length(mtime)]; % find gaps longer than a couple minutes
        chfirst = [1; chlast(1:end-1)+1];
        dchunk = chlast - chfirst;
        ii = dchunk >= 450;
        chlast  = chlast( ii);
        chfirst = chfirst(ii);
        
        for ichunk = 1:length(chfirst)
            iich = chfirst(ichunk):chlast(ichunk);
            if     any( isfinite(Urelsnd(:,iich)) )
                U = nanmean( Urelsnd(:,iich), 2);
%             elseif any( isfinite(Urelshp(:,iich)) )
%                 U = nanmean( Urelshp(:,iich), 2);
            else
                continue
            end
            epsilon(:, count+1) = dissipation_spectr_acov3( Dv(:,iich), U, fs );
%             pcolor(mtime(iich), height, Dv(:,iich)); shading flat; datetick('x'); colorbar; caxis(mean(nanstd(Dv(:,iich)))*[-3 3]); ylim([0, 1]);
%             hold on; semilogx(epsilon(:, count+1), height)
            Ulev(:, count+1) = U;
            time_start(count+1) = mtime( chfirst(ichunk) );
            time_end(  count+1) = mtime( chlast( ichunk) );
            count = count + 1;
        end
        
        %% plot cross sections by continuous chunk
        if plotit
            b2rcolormap(16);
            clf; ax = [];
            fig = gcf;
            ktop = 35;
            for ichunk = 1:length(chfirst)
                iich = chfirst(ichunk):chlast(ichunk);
                ax(ichunk) = subplot(length(chfirst),1, ichunk, 'align');
                t = mod(S.dechour(iich), 1)*60;
                pcolor(t, range(1:ktop), Dv(1:ktop,iich)); shading flat;
                caxis([-1.5 1.5]); colorbar();
                xlim(floor(t(1)/10)*10 + [0 10])
                xlabel('minute'); ylabel('height (km)')
                set(ax(ichunk), 'ylim', [0 1], 'fontsize',10)
                fig.Units = 'inches'; fig.OuterPosition = [1 1 5.5 7];
                if ichunk == 1
                    title(['Streamline XR lidar velocity anomalies ' datestr(mtime(chfirst(1)), 'yyyy-mm-dd HH')],...
                    'fontweight','normal')
                end
                % print('-dpng',['./plots/lidar_dopplervel_' flist(ifile).name(10:20) '.png']);
            end
        end
    end % hourly files
    fprintf(1, '\n')
end % daily folders

% truncate data to size
epsilon    = epsilon( :, 1:count );
Ulev       = Ulev(    :, 1:count );
time_start = time_start( 1:count );
time_end   = time_end(   1:count );

%save('dissipation', 'epsilon','time_start','time_end')
save('dissipation20201129', 'time_start','time_end','height', ...
    'epsilon', 'Ulev')

%% plots

% time gaps with NaN
t = paren([time_start time_end+2/1440]',':');
[eps, Ulv] = deal( NaN(size(epsilon).*[1 2]) );
eps(:, 1:2:(2*count)) = epsilon;
Ulv(:, 1:2:(2*count)) = Ulev;

b2rbrewcmap(12);
subplot(2,1,1)
pcolor(t, height, Ulv); shading flat;
set(gca,'color',0.7+[0 0 0], 'ylim',[0 1])
datetick('x')


b2rbrewcmap(12);

clf
subplot(2,1,1)

subplot(2,1,2)
pcolor(t, height, log10(eps)); shading flat;
hc=colorbar; hc.YLabel.String='log_{10}(\epsilon/m^2s^{-3})';
caxis([-6 -1]);
set(gca,'color',0.7+[0 0 0], 'ylim',[0 1],'fontsize',14)
xlim(datenum(2019,07,[17 18.6]))
xticks(datenum(2019,07,17:0.25:18.6))
datetick('x', 'HH', 'keeplimits','keepticks')
title('July 17-18', 'fontweight','norm')

% the horizontal winds go missing at the end of the record
subplot(2,1,1)
pcolor(H.time-datenum(2019,7,0), range, H.urel); shading flat
colorbar
caxis([-25 25])
subplot(2,1,2)
pcolor(H.time-datenum(2019,7,0), range, H.vrel); shading flat
colorbar
caxis([-25 25])

subplot(2,1,1)
pcolor(H.time-datenum(2019,7,0), height, H.uwnd); shading flat
colorbar
caxis([-25 25])
subplot(2,1,2)
pcolor(H.time-datenum(2019,7,0), range, H.vwnd); shading flat
colorbar
caxis([-25 25])

pcolor(time_start, height, log10(epsilon)); shading flat; datetick('x','keeplimits')
 
%{
% compare wind from lidar to from soundings
subplot(2,1,1)
pcolor(H.time-datenum(2019,7,0), height, H.uwnd); shading flat
colorbar
caxis([-25 25])
subplot(2,1,2)
pcolor(mtime-datenum(2019,7,0), height, usnd); shading flat
colorbar
caxis([-25 25])
%}
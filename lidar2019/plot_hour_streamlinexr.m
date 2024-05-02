% plot_hour_streamlinexr

nrange = 104; % number of range gates

%% pick a file to read
datadir = '/Volumes/MISOBOB/Halo Lidar sn0610-06 LEG 2/';
dd = dir(fullfile( datadir, '2019*' ));
dates = {dd.name};

% file = '/Volumes/MISOBOB/Halo Lidar sn0610-06 LEG 2/20190712/Stare_06_20190712_00.hpl';
id = 1;
hour = 0;
hh = sprintf('%02i', hour);
starefile = fullfile(datadir, dates{id}, ['Stare_06_' dates{id} '_' hh '.hpl']);
% scanfiles  = fullfile(datadir, dates{id}, ['User5_06_' dates{id} '_' hh '*.hpl']);

% read the stare file into the structure a
a = read_streamlinexr_stare(starefile, nrange);
% mask our where returns are too weak
thr = 1.5; % seems to work, by experimentation
mask = a.intensity > thr;
nanmask = zeros(size(mask));
nanmask(~mask) = NaN;

% subtract out mean heave
heave = nanmean( nanmask + a.dopplervel , 1 );

height = 30.0*(0:nrange-1)/1e3;
iih = height < 1.030;

ichend = [find(diff( 60 * a.dechour ) > 1); length(a.dechour)];
icstart = [1; ichend(1:end-1) + 1];

for ichunk = 1:6
clf; clear ax
ax(1) = subplot(2,1,1,'align');
imagesc(a.dechour*60, height(iih), log(a.intensity(iih,:)) ); axis xy
ax(2) = subplot(2,1,2,'align');
pcolor(a.dechour*60, height(iih), nanmask(iih,:) + a.dopplervel(iih,:) - heave ); axis xy
set(ax(:), 'ylim', [0 1],'fontsize',14)
xlabel('minute'); ylabel('height (km)')
orient landscape

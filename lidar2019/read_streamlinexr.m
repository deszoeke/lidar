function a = read_streamlinexr(file, nrange)

% file = '/Volumes/MISOBOB/Halo Lidar sn0610-06 LEG 2/20190712/Stare_06_20190712_00.hpl';

x = readtable(file, 'filetype','text', 'headerlines',17);

% get ray time and angles
% Decimal time (hours)  Azimuth (degrees)  Elevation (degrees) Pitch (degrees) Roll (degrees)
timeangles = x(1:(nrange+1):end,:);
a.dechour   = timeangles{:,1};
a.azimuth   = timeangles{:,2};
a.elevation = timeangles{:,3};
a.pitch     = timeangles{:,4};
a.roll      = timeangles{:,5};

nt = length(a.dechour);

% get ray data
idata = mod((1:size(x,1))'-1, nrange+1) > 0;
% Range Gate  Doppler (m/s)  Intensity (SNR + 1)  Beta (m-1 sr-1)
a.rangegate  = reshape(x{idata,1}, [nrange, nt]);
a.dopplervel = reshape(x{idata,2}, [nrange, nt]);
a.intensity  = reshape(x{idata,3}, [nrange, nt]);
a.beta       = reshape(x{idata,4}, [nrange, nt]); % backscatter cross section [scatter area / sample volume / steradian]



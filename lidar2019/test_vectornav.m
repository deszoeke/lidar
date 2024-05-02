% test_vectornav.m

switch computer
    case 'PCWIN64'
        addpath(genpath('C:\Users\Simon de Szoeke\Documents\MATLAB\user_tools\graphics\'));
    case 'MACI64'
        addpath(genpath('/Users/sdeszoek/Documents/MATLAB/user_tools/graphics/'));
end

vndir = '~/Data/cruises/MISOBOB_2019/SR1911/lidar/HaloLidar_sn0610-06_leg2/Stable_Data/';

flist = dir(fullfile(vndir, 'VectorNavData*'));

anom = @(x) x - nanmean(x);

% from shell> wc VectorNavData* --> 55896040 total
nt =  55896040; % about a month at 20 Hz
dnum = NaN(nt,1); % allocate
ct = 0;
for ifile = 1:length(flist)
    flist(ifile).name
    mdate(ifile) = datenum(flist(ifile).name(14:32), 'yyyy_mm_dd_HH_MM_SS');
    
    % awfully slow
%     T = readtable( fullfile(vndir, flist(ifile).name), 'filetype','text', 'delimiter',',' );
    % load is also slow, and eventually crashes on parsing formatted text after
    
    fid = fopen( fullfile(vndir, flist(ifile).name) );
    fr = fileread( fullfile(vndir, flist(ifile).name) ); % char array, slow
    i=0;
%         lin = fgetl( fid );
%         a = textscan( fid ,['%*s' repmat('%f', 1,20)], 'delimiter',',' );
%         strlen = length( a{1} );
%         b = textscan( fid ,['%*s' repmat('%f', 1,20)], 'delimiter',',', 'collectoutput',1 ); % slow
%         T = cat(1, cell2mat(a(2:end)), b{:}); % also slow
%         t0 = datenum( T{1,1}{:}(1:end-6), 'ddd mmm dd HH:MM:SS yyyy' );
    T = textscan( fid ,['%*s' repmat('%f', 1,20)], 'delimiter',',', 'collectoutput',1 );
    
    % matlab UTC date number from GPS time in ns.
    % 18 s offset is according to GPS standard
    UTC = datenum(1980, 1, 6, 0, 0, T(:,1)*1e-9 + (37-19) );
    % linear vel from integrating accelerations
    vel = cumsum(T(:,15:17)-mean(T(:,15:17)))*0.05; % component 3 is vertical
    % note velocities start in columns 18:20 when GPS initializes
        
%     for it = 1:nc
%         dnum(ct+it) = datenum( T{it,1}{:}(1:end-6), 'ddd mmm dd HH:MM:SS yyyy' );
%     end
    ct = ct + nc;
end

%{
columns of T
 1      "GpsTime", nanoseconds since 1980-01-06 00:00:18
 2: 4   "Yaw", "Pitch", "Roll",
 5: 8   "Quat0", "Quat1", "Quat2", "Quat3",
 9:11   "Latitude", "Longitude", "Altitude", % zeroed
12:14   "MagNED0", "MagNED1", "MagNED2",
15:17    "LinAcc0", "LinAcc1", "LinAcc2",
18:20    "VelNED0", "VelNED1", "VelNED2");

    % The GPS time jumps when the linear velocities turn on.
    % There?s a jump of 1.25*10^18 ns = 39.5 years in the 17-Jun-2019 06:17:11 file.
    % The GPS must assumes it?s still 1980 until then.

    % first time stamp matches filename.
    % -0400 maybe corresponds to Eastern daylight time zone (GMT-4h)
    % used in South Bend, Indiana. So time is in UTC, and the -0400 means
    % subtract 4 to get Eastern time.
    
    % frequency is 0.05 s, 20 Hz
    plot(T{:,3:9},'.')
    plot(T{:,13:19},'.')
    % table T(158,21) variables
    %  2     is monotonic, time in nanoseconds, 10^-9 s, dt = 0.05 s
    %  3: 5  are O(10) to ~100
    %  6: 9  are O(0.1),
    % 10:12  are exactly zero
    % 13:15  are a little bit noisy
    % 16:18  are noisiest
    % 19:21  are zero

VectorNav::VectorNav() : vs() {
  shared_memory = new VNSM();
  shared_memory->name(
    "GpsTime",
    "Yaw", "Pitch", "Roll",
    "Quat0", "Quat1", "Quat2", "Quat3",
    "Latitude", "Longitude", "Altitude",
    "MagNED0", "MagNED1", "MagNED2",
    "LinAcc0", "LinAcc1", "LinAcc2",
    "VelNED0", "VelNED1", "VelNED2");
data_logger<uint64_t, float, float, float, float, float, float, double, double, double, float, float, float, float, float, float, float, float, float, float, float> vn_logger(vnsm, vn_data_path, 20_hz);


%}

% 1. add heave to Doppler velocity
% 2. (optional) regress out horizontal component times cos,sin tilt velocity
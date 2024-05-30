function [ mld, kmld0, kmld1, epsbar, epstop, mldgrad, kgrad ] = calc_mld( epsilon, epy, zData, xrthr )
% [ mld, kmld0, kmld1, epsbar, epstop, mldgrad, kgrad ] = calc_mld( epsilon, epy, zData, xrthr )
% mld: ML depth from epsilon located in the quiescent air above the turb. ML
% 1. Calculate mixed layer depth mld from the filtered dimensional quantity epy (linear not dB),
% and height variable zData, as the level where epy drops to 1/xrthr of its
% vertical mean below and including that level.
% 2. Find the max gradient
% 3. Refine the max gradient by interpolating the inflection point.
% Deal with missing points properly!
%  Simon de Szoeke :: 2020 :: HRDL dissipation calculation

%% 1. Find where dissipation first is less than 1/3 mean dissipation below
% that level.
% Also find kmld, first index where this is true.
% kmld0 is the first level where dissipation < 1/3 mean
% below.
[nt, kmax] = size( epy );
mld    = NaN(nt, 1); % 4040 turbulence profiles
epsbar = NaN(nt, 1);

% epy is the filled and filtered dissipation
iiy = isfinite(epy);
meandiss = cumsum(epy,2)./cumsum(iiy,2);
stddiss  = cumsum((epy-meandiss).^2)./cumsum(iiy,2).^(3/2); % std of the mean
xr = meandiss ./ epy;

for it = 1:nt
    strev = find(xr(it, end:-1:1) > xrthr, 1,'last'); % search downward for "start" index
    st = kmax - strev;
    strange = (max(1,st-1)) : min(kmax,st+1);
    xrf = xr(it,strange);
    strange = strange(isfinite(xrf));
    if length(strange) > 1
        mld(it)    = interp1( xr(it,strange), zData(strange), xrthr  , 'linear' );
        epsbar(it) = interp1( zData(strange), meandiss(it,strange), mld(it), 'linear' );
    end
end

% kmld0 is the first level where dissipation < 1/3 mean
kmld0 = floor( (mld-zData(1)) / (zData(2)-zData(1)) ) + 1;
% regardless of whether there is a real valid dissipation estimate there

epstop = zeros(size(kmld0,1), 1);
ii = kmld0 > 0;
epstop(ii) = diag( epy(ii,kmld0(ii)) );


%% 1.5 Move kmld0-->kmld1 _up_ to first estimated dissipation, so not to locate MLD at
% an interpolated height.
kmld1 = kmld0;
shift = zeros*kmld0;
iishift = diag(epsilon(kmld1>0, kmld0(kmld1>0)))<=0;
for it = find( iishift )
    shift(it) = find( epsilon(it,kmld0(it):kmax) > 0, 1, 'first') - 1;
end
kmld1 = kmld0 + shift;
mld0(shift>0) = zData(kmld0(shift>0)); % also move mld0 to valid epsilon


%% 2. get kgrad and initial mldgrad0 = zData(kgrad) where gradient is minimum (neg.)
% zoff = zData(1:kmax)' - mld0;
% epstop = zeros(size(zoff,1),1);
ii = isfinite(kmld1);
epstop = NaN(size(kmld1));
epstop(ii) = diag( epy(ii, kmld1(ii)) );
% kgrad < kmld1
iinml = 1:kmax <= kmld1 + 1;              % is below k of MLD top
iimin = islocalmin( diff(epy, 1, 2), 2 ); % is local min of gradient of filtered dissipation
% find the highest min vertical derivative of epy below kmld for each time
% !Searching for highest minimum derivative is probably not very stable! %
kgrad = zeros(length(kmld0),1, 'int8');
for it = 1:length(kgrad)
    tmp = find( iinml(it, 1:end-1) & iimin(it, :), 1, 'last' );
    if ~isempty(tmp)
        kgrad( it ) = tmp;
    end
end
mldgrad0 = zeros(size(kgrad));
dz = zData(2) - zData(1);
jj = isfinite(kgrad) & (kgrad>0);
mldgrad0(jj) = zData(kgrad(jj)) + dz/2; % 0th guess picks maximum (neg) gradient layer

%% 3. refine the gradient MLD by interpolating the inflection point
% find max gradient in dissipation by interpolating inflection point
mldgrad = mldgrad0;
for it = 1:length(kgrad)
    % stencil in neighborhood of kgrad (was kmld)
    stencil = max(1, kgrad(it) - 4) : min(kgrad(it) + 1, kmax);
    % exclude interpolated parts when epsilon not retrieved
    iisf = epsilon(it,stencil) > 0;
    stencil = stencil( stencil>=1 & stencil<=kmax & iisf )';
    % interpolate gradients between missing epsilon
    dz1 = diff( zData(stencil) );
    zm = zData(stencil(1:end-1)) + dz1;
    de1 = diff( epy(it, stencil ), 1, 2 )' ./ dz1;
    %2d order derivative along 2d dimension
    de2 = diff( de1 ) ./ diff( zm );
    ii2 = isfinite( de2 );
    iist = stencil(2:end-1);
    if length(ii2) >= 2
        mldgrad(it) = interp1( de2(ii2), zData( iist(ii2) ), 0 );
    end
end

%{

% To compare the dissipation at k to the second highest dissipation below
% k, store the second highest value in a register.

% Running register reg keeps length(reg) largest values in ascending order
paren = @(x, varargin) x(varargin{:}); % for indexing or splatting arguments into functions
reg = zeros(2,1);
update_max_reg = @(reg, x) paren( sort([x; reg]), 2:length(reg)+1 );

%}

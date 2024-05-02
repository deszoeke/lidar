cd('/Users/sdeszoek/Data/cruises/MISOBOB_2019/SR1911/scientific_analysis/programs/lidar');

addpath('/Users/sdeszoek/Documents/MATLAB/user_tools/thermo/');
addpath('/Users/sdeszoek/Data/cruises/MISOBOB_2019/SR1911/scientific_analysis/programs/flux');
% flux data
[f10, varname, description] = get_flux_10min();
% lidar TKE dissipation data
load('dissipation20201129', 'time_start','time_end','height', ...
    'epsilon', 'Ulev');

sfc_eps = (0.52 * f10.cu).^(3/2);
[bflx, bflx_s] = buoyancy_flux(f10.hsb, f10.hlb, f10.ta, f10.qa, f10.Pmb);

% time gaps with NaN
t = paren([time_start time_end+2/1440]',':');
[eps, Ulv] = deal( NaN(size(epsilon).*[1 2]) );
eps(:, 1:2:(2*size(epsilon,2))) = epsilon;
Ulv(:, 1:2:(2*size(epsilon,2))) = Ulev;

%% plots
b2rbrewcmap(12);

clf
ax(1) = subplot(3,1,1);
area(datenum(2019,1,f10.doy), f10.rs/1e2, 'edgecolor',0.5*[.9 .8 .2],'facecolor',[.9 .8 .2]); % mm/h
hold on
area(datenum(2019,1,f10.doy), f10.rain/10, 'edgecolor','none','facecolor',0.8*[.7 1 .7]); % mm/h
hold on
plot(datenum(2019,1,f10.doy), f10.wspd_scs, 'k', 'linewidth',1.4)
yticks(ax(1), 0:2:10); %yticklabels([0 5 10 25 30]);
xlim(datenum(2019,07,[18 23]))
xticks(datenum(2019,07,18:23))
datetick('x', 'dd', 'keeplimits','keepticks')
cb(1) = colorbar(ax(1)); cb(1).Visible = 'off';
ylabel(ax(1), {'rain (10 mm/h)', 'wind (m/s)'})

ax(2) = subplot(3,1,2);
plot(datenum(2019,1,f10.doy), f10.ta, 'b','linewidth',1.4)
hold on
plot(datenum(2019,1,f10.doy), f10.ts, 'r','linewidth',1.4)
xlim(datenum(2019,07,[18 23]))
xticks(datenum(2019,07,18:23))
datetick('x', 'dd', 'keeplimits','keepticks')
cb(2) = colorbar(ax(2)); cb(2).Visible = 'off';
ylabel(ax(2), 'T (°C)')

ax(3) = subplot(3,1,3);
plot(datenum(2019,1,f10.doy), [bflx bflx_s], 'linewidth',1.4)
hold on
plot(datenum(2019,1,f10.doy([1 end])), zeros(2,1),'k')
xlim(datenum(2019,07,[18 23]))
xticks(datenum(2019,07,18:23))
datetick('x', 'dd', 'keeplimits','keepticks')
cb(3) = colorbar(ax(3)); cb(3).Visible = 'off';
ylabel(ax(3), 'bouyancy flux (m^2 s^{-3})')

set(ax(:),'fontsize',14)
set(ax(:),'xlim', datenum(2019,07,[18 23]))

set(gcf,'inverthardcopy','off','color','w')
% print -depsc2 Jul18_23met.eps


%% add TKE dissipation

plot(ax(1), datenum(2019,1,f10.doy), [bflx bflx_s]*1e4, 'linewidth',1.4)

ax(3) = subplot(3,1,3);
hold off
pcolor(t, height, log10(eps)); shading flat;
hc=colorbar; hc.YLabel.String='log_{10}(\epsilon/m^2s^{-3})';
caxis([-6 -3]);
set(gca,'color',0.7+[0 0 0], 'ylim',[0 1],'fontsize',14)
xlim(datenum(2019,07,[18 23]))
xticks(datenum(2019,07,18:23))
datetick('x', 'dd', 'keeplimits','keepticks')
ylabel(ax(3), 'height (km)')

% print -depsc2 Jul18_23met_tkediss.eps
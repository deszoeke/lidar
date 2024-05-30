% Calculate cloud top TKE dissipation from a vertical structure function
% for VOCALS 2008 W-band radar
% Simon de Szoeke 2013.04.26
% 2018.07

%% set environment and constants
%cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband
%cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/wband
warning('off','MATLAB:interp1:NaNinY')
run('../../read_parameters')
machine=char(java.net.InetAddress.getLocalHost.getHostName);
if regexp(machine,'squall')==1
    stem='~/Data/cruises/VOCALS_2008/RHB';
    way_Pinsky_retrieval='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/wband/Pinsky/retrieval/';
elseif regexp(machine,'fog')==1
    stem='/Volumes/Pegasus2/cruises/VOCALS_2008/RHB';
    way_Pinsky_retrieval='/Volumes/Pegasus2/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky/retrieval/';
end

yyyy=sprintf('%04d',year);
stday=318;
enday=335;
fs=3.5; %Hz
crosscomp=3/4;
kolmogorov=0.54; % Matches atmospheric boundary layer estimates and Sreenivasan 1995
% kolmogorov=0.5; % probably only 1 digit of precision, Sreenivasan 1995
C1prime=4/3*kolmogorov; % as in Pope eqn 6.243
factr=crosscomp/kolmogorov; % 1/C1prime, used for my dissipation calculation

% universal constant for 2nd order structure function
% longitudinal structure function constant C2 from Pope 6.2 (p.193) after Saddoughi and Veeravalli (1994)
C2ll=2.0; % 2.1, up to 2.2? Wiles et al. 2006
factrz=1/C2ll;
factrx=3/4/C2ll;
% Kolmogorov (viscous) scale for atmospheric boundary layer is
% eta=(nu^3/epsilon); epsilon~3e-4 m^2/s^2/s, kinematic viscosity nu=2e-5 m^2/s
% --> eta= 2.3 mm
%     25 m = 1000 eta;  1000 m = 440 000 eta

dten=10*60; % seconds per window (10 min)
diten=floor(fs*dten); % samples per window

%% structure function parameters
nmin=100; % min number of Doppler velocities to calc structure function in a 10 min interval

%% load 10-min ship-relative wind at cloud top interpolated from soundings
if regexp(machine,'squall')==1
    stem='~/Data/cruises/VOCALS_2008/RHB';
    load([stem '/Scientific_analysis/programs/wband/uv_cloudtop10.mat']);
    load([stem '/Scientific_analysis/programs/wband/uv10radar.mat']);
elseif regexp(machine,'fog')==1
    stem='/Volumes/Pegasus2/cruises/VOCALS_2008/RHB';
    load([stem '/Scientific_analysis/programs/VOCALS2008_programs/wband/uv_cloudtop10.mat']);
    load([stem '/Scientific_analysis/programs/VOCALS2008_programs/wband/uv10radar.mat']);
end
% ship maneuvers removed
% horizontal ship-relative speed
% U=sqrt(utoprel.^2+vtoprel.^2); % m/s interpolated relative U from soundings and ship velocity
% Uradar=sqrt(uradarrel.^2+vradarrel.^2);

%% loop length numbers
nhr=24*(enday-stday+1);
n10min=nhr*6;
%% allocate output data variables
jmax=83; lmax=jmax-1; % 72,71 % actually jmax could be jmax-1 too
ydayw10=NaN+zeros(n10min,1);
ktop=zeros(n10min,1,'int8');
kmax=zeros(n10min,1,'int8');
count=zeros(n10min,jmax,lmax,'int16');
D=NaN(n10min,jmax,lmax);
[noise,A]=deal(NaN(n10min,lmax));

%% loop over files by day, hour
for iday=stday:enday
    fprintf(1,'\n%3d ',iday)
    for idhr=0:23
        fprintf(1,' %02d ',idhr)
        retrfile=dir([way_Pinsky_retrieval sprintf('Pinsky_wRetrieval%04i_%03i_%02i.mat',2008,iday,idhr)]);
        if ~isempty(retrfile)
            load([way_Pinsky_retrieval retrfile.name]);
            nt=size(W,1);
            % calculate cloud-top dissipation from inertial range of Pinsky-adjusted vel spectra
            % loop through, interp and calc spectra for 10-min chunks
            for ind=1:round(nt/diten) % index of the 10 minute interval in the hour
                bigind=24*6*(iday-stday)+6*idhr+ind;
                i10=(ind-1)*diten+1; % first index in the 10 minute interval
                ydayw10(bigind,1)=iday+(idhr+(ind-1)/6)/24; % current time of this 10-minute window
                %Utop=double(interp1(yday10,U,ydayw10(bigind)+5/(24*60))); % interp ship-rel wind speed to the middle of the current time window
                % find the max coverage of nonnoise returns --> nmax
                ii10=i10:min(i10+diten-1,nt);
                nw=sum(isfinite(W(ii10,:)));
                [nmax,knmax]=max(nw);
                % choose the highest vertical level ik that has >0.5*nmax of this coverage
                ik=find(nw>100,1,'last');  % highest
                hk=find(nw>100,1,'first'); % lowest
                    
                % compute dissipation by longitudinal 2nd order stucture function
                % if there is 10 min cloud Doppler velocity data (A cloud top is identified.)
                if ~isempty(ik) && ik<100 % A cloud top is identified.
                    ktop(bigind)=min(ik+1,100);
                    ktop(bigind)=min(ik+1,100);
                    kk=(max(hk-1,4):ktop(bigind))';
                    w10=W(ii10,kk); % subset w10 and filter out extreme values
                    w10(bsxfun(@gt,abs(bsxfun(@plus,w10,-nanmean(w10))),4*nanstd(w10)))=NaN;
                    %iiwt=isfinite(w10);
                    nj=length(kk);
                    kmax(bigind)=nj;
                    for j=1:nj % could be 1:nj-1 % index relative to 10min window cloud top k=ktop+(j-1)
                        % k1=nj-j+1; % indexes w10( ,k1)
                        iij0=isfinite(w10(:,nj-j+1)); % j=1 -> cloud top
                        for l=1:nj-j % structure fcn displacement index
                            iijl=isfinite(w10(:,nj-j+1-l));
                            iijj=iij0 & iijl;
                            count(bigind,j,l)=sum(iijj);
                            dw=w10(iijj, nj-j+1)-w10(iijj, nj-j+1-l);
                            D(bigind,j,l)=mean(dw.*dw);
                        end
                    end
                end % there is 10 min Doppler velocity data                
            end % 10 min
        end % there is hourly data
    end % hr
end % day

% truncate data
ydayw10=ydayw10(1:bigind);
ktop=ktop(1:bigind);
kmax=kmax(1:bigind);
jmax=max(kmax); lmax=jmax-1;
count=count(1:bigind,1:jmax,1:lmax);
D=D(1:bigind,1:jmax,1:lmax); % D(time, height, displacement)

%% calculate dissipation from the longtudinal structure function
% preliminaries
dh1=double(h(50))-double(h(49));
m=bsxfun(@plus,(0:jmax-1)',1:lmax); % m=j-1+l; midpoint index
mid=double(m)*0.5*dh1;                % midpoint displacement below cloud top (m)
r23=(dh1*(1:double(lmax))').^(2/3);
r23_2=repmat(shiftdim(r23,-1),[jmax 1]); % 2d version of r^(2/3)
% solve D = N + A*r^(2/3) for noise N and coefficient A
% x=lsqlin(Cee,Dee,Aye,bee) minimizes the residual Cee*x-Dee subject to A*x<=b
Aye=-eye(2); % limiting conditions for solution
bee=zeros(2,1);

fii=find(sum(sum(isfinite(D),2),3)>100); % times that have reasonable ensembles of structure functions

options = optimoptions('lsqlin','Algorithm','interior-point');
for itime=1:length(fii) % pick a time
    nj=kmax(fii(itime));
    D2=squeeze(D(fii(itime),:,:));
    count2=squeeze(count(fii(itime),:,:));
%     jjll=isfinite(D2);
    jjll=count2>10;
    for pickm=1:nj-1 % compute structure function for same midpoint height
        jjllmm = jjll & m==pickm;
        if sum(jjllmm(:))>5
            Dee=D2(jjllmm); % vectorize
            Cee=[ones(sum(jjllmm(:)),1) r23_2(jjllmm)];
            x=lsqlin(Cee,Dee,Aye,bee); % x is [noise A]
            noise(fii(itime),pickm)=x(1);
            A(fii(itime),pickm)=x(2);
        end
    end
end
% x(2)=0 probably is bad
% lots of identically valued A and identically valued noise.

epsz=(A*factrz).^1.5;
semilogy(ydayw10,epsz,'.')

pcolor(ydayw10,-mid(1,:),log10(epsz)'); shading flat;
caxis([-5 -2.5])
set(gca,'color',0.7+[0 0 0],'fontsize',16)
colorbar
title('dissipation from longitudinal structure function log_{10}[m^2 s^{-3}]','fontweight','normal')
xlabel('yearday 2008')
ylabel('midpoint displacement from cloud top [m]')
saveas(gcf,'epsz_cloudtoprel.eps','epsc')

% reconstruct dissipation on absolute height grid
dh2=dh1/2;
habs=flipud(h(max(ktop))-dh2*(1:size(epsilonz,2))'); % increasing
epsilonz=NaN(length(ydayw10),max(ktop)+size(mid,2)-1);
fl=find(isfinite(epsz)); % going in the direction that sparse is better than full, but not quite there yet
for fli=1:length(fl) % very fast
    l=fl(fli);
    [i,em]=ind2sub(size(epsz),l);
    epsilonz(i,double(ktop(i))-em) = epsz(l);   
end

b2rcolormap(17)
pcolor(ydayw10,habs,log10(epsilonz)'); shading flat;
caxis([-5 -2.5]); ylim([400 1300])
colorbar; set(gca,'color',0.7+[0 0 0])

% compare with spectral dissipation calculation
% at every other vertical level of the strfcn calc.

load('radarturb3.mat')
find(abs(h-habs(end-1))<10)
%habs(125)==h(76);
% koffs=(length(habs)-1) - find(abs(h-habs(end-1))<10); % offset at top
koffs=find(abs(h-habs(1))<5); % offset at base, 1-based, =14
% habs(1) == h(koffs)
plot(koffs-1+(1:length(habs)/2),habs(1:2:end),'.')
nk=fix(length(habs)/2); % 63
% sum( abs( h(koffs+(0:nk-1)) - habs(1:2:end) )>5 ) % 0
% h(koffs+(0:nk-1))==habs(1:2:end)
ii=(any(isfinite(epsilonz),2));
x=epsilonk(ii,koffs+(0:nk-1));
y=epsilonz(ii,1:2:end);
loglog(x(:),y(:),'.','markersize',1)
hold on
plot([1e-6 1],[1e-6 1],'k-')
plot([1e-6 1],2*[1e-6 1],'k--')
ylabel('structure function dissipation')
xlabel('spectral dissipation')
set(gca,'fontsize',14)
axis([1e-6 1e-2 1e-6 1e-2])
%[nanmedian(x(:)) nanmean(x(:)); nanmedian(y(:)) nanmean(y(:))]
% not well correlated, but distributions agree pretty well in the
% mean,median
% strfcn dissipation is not well localized, and I haven't conditioned on
% regions of nearly constant dissipation.

% ...or could do by interpolation
% xi=find(any(isfinite(epsilonz),2));
% interp2(xi,habs


plot([0; 120],x(1)+x(2)*[0; 120],'r-')
hold on
plot(r23_2(jjllmm),Dee,'.')

clf
pcolor(dh1*(1:60),dh1*(0:50-1),squeeze(D(fii(itime),1:50,1:60))); shading flat;
cx=caxis();
hold on
contour(dh1*(1:60),dh1*(0:50-1),mid(1:50,1:60),0:100:1000,'k')
set(gca,'tickdir','out','color',0.7+[0 0 0],'fontsize',14)
xlabel('displacement (m)'); ylabel('high point distance below cloud top (m)')
caxis(cx);
colorbar('East')

% QC for enough points in the spectrum
isenuf=count>1800;
isenufk=countk>1800;

% Already QCd Pinsky retrieval for missing motion data and motion
% correction failure, (85% good) and restored some motion values using
% the Crossbow accelerometer/gyro.
%  and, a posteriori, variance <0.4 m^2/s^2
isgdtmp=isenuf;
% nanstd(sigv2(isgdtmp))*3 ~ 0.4
isgd=isgdtmp & sigv2<0.4;
sp=spec(isgd,:);
ep=epsilon(isgd);
yd=ydayw10(isgd);

spk=speck(isgd,:);
epk=epsilonk(isgd,:);


set(gca,'ydir','normal')
%}

%% save data
save cloudtopturb6.mat ydayw10 F spec epsilon sigv2 count knoisetop Snoisetop Sthreshtop
% vertically-resolved data
save radarturb6.mat ydayw10 F speck epsilonk sigv2k countk knoisek Snoisek Sthreshk

%% plot nondimensional spectra scaled by Kolmogorov power spectral density (dissipation*dynvisc^5)^1/4
% for cloud top
dynvisc=1.63e-5; %m^2/s, at ~900 hPa
Kpsd=sqrt(sqrt(epsilon*(dynvisc^5)));  % Kolmogorov power spectral density (dissipation*dynvisc^5)^1/4  [U^2 / L^-1]
Kks=sqrt(sqrt(epsilon./(dynvisc^3)));  % Kolmogorov wavenumber (dissipation*dynvisc^3)^1/4  [L^-1]
eta=1./Kks;                            % Kolmogorov length scale
wavenumber=2*pi*bsxfun(@rdivide,F',utop);
kscaled=bsxfun(@rdivide,wavenumber,Kks); % nondimensional wavenumber is scaled by Kolmogorov wavenumber
speck=bsxfun(@times,spec,utop/(2*pi));   % wavenumber spectrum (spec is frequency spectrum!)
Sfscaled=bsxfun(@rdivide,spec,Kpsd);           % nondimensional (TKE) power frequency-spectral density scaled by Kolm'v PSD
Skscaled=bsxfun(@times,Sfscaled,utop/(2*pi));  % nondimensional (TKE) power wavenumber-spectral density scaled by Kolm'v PSD

% cut off noise
Scut=Skscaled;
for i=1:length(knoisetop)
    Scut(i,max(1,knoisetop(i)):end)=NaN;
end
Scutadj=bsxfun(@minus,Scut,Snoisetop.*utop/(2*pi)); % subtract noise
Sadj=bsxfun(@minus,Skscaled,Snoisetop.*utop/(2*pi)); % subtract noise

% plot scaled not by energy containing scale or by noise, but by unseen inferred Kolmogorov scale
clf
loglog(kscaled',Skscaled','.','markersize',1,'color',0.7*[1 1 1])
hold on
loglog(kscaled(:,4:end)',Scut(:,4:end)','b.','markersize',1)
% loglog(kscaled(:,4:end)',Scutadj(:,4:end)','c.','markersize',1)
plot([5e-5 1e-2],C1prime*[5e-5 1e-2].^(-5/3),'r-','linewidth',1.4)
set(gca,'fontsize',16,'xlim',[1e-5-10*eps 1e-1])
title('Kolmogorov-scaled w wavenumber spectra','fontweight','normal')
ylabel('power spectral density   S_{ww}(k)/(\epsilon\nu^5)^{1/4}')
xlabel('wavenumber   k\eta')

% plot compensated spectrum
% clf
loglog(kscaled',Skscaled'.*kscaled.^(5/3)','.','markersize',1,'color',0.7*[1 1 1])
hold on
loglog(kscaled(:,4:end)',Scut(:,4:end)'.*kscaled(:,4:end).^(5/3)','b.','markersize',1)
plot([1e-5 1e-1],C1prime*[1 1],'r-')
%%
% saveas(gcf,'Kolmogorov_spectrum.eps','epsc')

errorbar(nanmean(Scut(:,4:end)'.*kscaled(:,4:end).^(5/3)',2),nanstd(Scut(:,4:end)'.*kscaled(:,4:end)'.^(5/3),2)./sqrt(sum(isfinite(Scut(:,4:end)'.*kscaled(:,4:end)'),2)))
% k,kscaled changes with the wind, so binavg
kbin=logspace(5e-5,5e-2,31);
[mb,sb,wb,nb]=binavg(kbin,kscaled(:,4:end),Scut(:,4:end).*kscaled(:,4:end).^(5/3));
[mlb,slb,wlb,nlb]=binavg(kbin,kscaled(:,4:end),log(Scut(:,4:end).*kscaled(:,4:end).^(5/3)));
errorbar(kbin,mb,sb./sqrt(nb))
set(gca,'xscale','log')

%%% plot

%% test QC procedure
plot(ydayw10(isgd),sigv2(isgd),'r.')
hold on
plot(ydayw10,sigv2,'b')
plot(ydayw10,epsilon/1e-2,'k')
plot(ydayw10(isgd),epsilon(isgd)/1e-2,'m.')
set(gca,'ylim',[0 0.35])

%% scatter plot epsilon and power
clf
plot(epsilon(isgd),sigv2(isgd),'.')
hold on
set(gca,'yscale','log','xscale','log')
axis([0 0.01 0 1])
plot([1.7e-3 1.7e-2].^1.5,[2e-2 2e-1],'r') % epsilon=power^1.5 dependence
set(gca,'fontsize',14)
xlabel('dissipation')
ylabel('<w''w''>')

%% TKE dissipation time scale (estimated from spectra
% maybe long eddy-containing eddies truncated
clf
plot(ydayw10,sigv2./epsilon/60)
hold on
plot(ydayw10(isgd),sigv2(isgd)./epsilon(isgd)/60,'m.')
set(gca,'fontsize',14)
title('TKE dissipation timescale')
ylabel('minutes')
xlabel('2008 yearday')
% compare to PBL cloud top

%% time series
subplot(2,1,1)
x=ydayw10(isfinite(ydayw10)); y=epsilon(isfinite(ydayw10));
isx=isgdtmp(isfinite(ydayw10)) & sigv2(isfinite(ydayw10))<0.4;
% y(~isx)=NaN;
% plot(x-5/24,y)
y(~isx)=-.01e-3;
area(x-5/24,y,'facecolor',[0.5 0.6 0.7],'edgecolor','none')
hold on
plot(ydayw10(isgd)-5/24,epsilon(isgd),'.','markersize',4)
plot(ydayw10-5/24,5e-4*max(0,-cos(2*pi*(ydayw10-5/24))))
x6=x(1:6:end);
isxx=[isx; zeros(6-mod(length(isx),6),1)]; % pad end
y6=all(~isxx([1:6:end; 2:6:end; 3:6:end; 4:6:end; 5:6:end; 6:6:end]))';
plot(x6-5/24,.7e-3*y6-.1e-3,'r.')
set(gca,'tickdir','out','ylim',[0 1.5e-3],'fontsize',14,'xminortick','on')
xlabel('2008 yearday')
ylabel('TKE dissipation m^2 s^{-3}')
% saveas(gcf,'../CTdissipation/dissip_timeseries.eps','epsc')

%% QC outcome displayed diurnally
plot(mod(ydayw10*24-5,24),epsilon,'rx') % all points
hold on
plot(mod(ydayw10(isgd)*24-5,24),epsilon(isgd),'.') % only good points

% histogram
subplot(2,1,1)
hist(log10(epsilon),-5.5:0.05:-0.5)
title('all')
set(gca,'xlim',[-5.5 -0.5])
subplot(2,1,2)
hist(log10(ep),-5.5:0.05:-0.5)
set(gca,'xlim',[-5.5 -0.5])
title('QCd')
xlabel('log_{10}(\epsilon)')
% dissipation is lognormally distributed
% QCing based on count and motion correction has little effect on its distribution

%% shade spectra and plot 10x<w'w'>
ii=isfinite(sigv2);
b2rcolormap(15)
pcolor(ydayw10(ii),F,log10(spec(ii,:)')); shading flat
hold on
plot(ydayw10,10*sigv2,'k')

%% sort the spectra by power by power in their low, high freq ends
clf
subplot(2,1,1)
[ssp,iord]=sort(sp(:,2));
pcolor(log10(sp(iord,:))'); shading flat
hold on
plot(1e4*ep(iord),'k')
set(gca,'xlim',[1 length(iord)],'ylim',[0 33])

subplot(2,1,2)
[ssp,iord]=sort(sp(:,end-1));
pcolor(log10(sp(iord,:))'); shading flat
hold on
plot(1e4*ep(iord),'k')
set(gca,'xlim',[1 length(iord)],'ylim',[0 33])

%% divide spectra into modes with an EOF analysis
ii=isfinite(sigv2);
X=log10(spec(ii,2:end-1))-repmat(mean(log10(spec(ii,2:end-1))),[sum(ii) 1]); % de-mean log10 power series
[u,s,v]=svd(X);

subplot(2,2,1)
pcolor(X'); shading flat
title('full matrix')
subplot(2,2,2)
p=zeros(size(s)); p(1,1)=s(1,1);
pcolor((u*p*v')'); shading flat
title('mode 1')
hold on
plot(v(:,1)/max(v(:,1))*1e3,1:nwin/2-2,'color','k','linewidth',1.4)
subplot(2,2,3)
p=zeros(size(s)); p(2,2)=s(2,2);
pcolor((u*p*v')'); shading flat
title('mode 2')
hold on
plot(v(:,2)/max(abs(v(:,2)))*5e2+5e2,1:nwin/2-2,'color','k','linewidth',1.4)
subplot(2,2,4)
p=zeros(size(s)); p(3,3)=s(3,3);
pcolor((u*p*v')'); shading flat
title('mode 3')
hold on
plot(v(:,3)/max(abs(v(:,3)))*5e2+5e2,1:nwin/2-2,'color','k','linewidth',1.4)

%% composite diurnal cycle of dissipation on local hour: mean, median histogram
indx=bin2(0:24,mod(yd*24-5,24)); % local hour index 1:24; 1 corresp to 0-1 local
ledges=-5.5:0.1:-0.5;
edges=0:1e-4:1e-2;
for ih=1:24
    ii=indx==ih;
    epse(ih)=nanmean(ep(ii));
    epsm(ih,:)=quantile(ep(ii),[.25 .50 .75]);
    % histogram
    hlep(ih,:)=histc(log10(ep(ii)),ledges);
    hep(ih,:)=histc(ep(ii),edges);
    % histogram with median diurnal cycle removed
    y=ep(ii)-epsm(ih,2);
    hdlep(ih,:)=histc(log10(y),ledges);
    hdep(ih,:)=histc(y,edges-2e-3);
end

clf
plot(mod(yd*24-5,24),ep*1e3,'.','color',[0.7 0.7 0.7])
hold on
plot(0.5:23.5,epse*1e3,'ko') % local hour
lh=plot(-0.5:24.5,epsm([end 1:end 1],:)*1e3,'k'); % local hour
set(lh(2),'linewidth',2)
set(gca,'fontsize',16,'xlim',[0 24],'xtick',0:6:24,'xticklabel',[0:6:18 0])
ylabel('TKE dissipation (10^{-3} m^2 s^{-3})')
xlabel('local hour')
set(gca,'ylim',[0 1.5])
%print -depsc CTdissipation/dissipation_diel.eps

colormap(1-gray(16))
% log scale diurnal histogram
subplot(2,1,1)
pcolor(0:36,ledges,hlep([1:end 1:end/2+1],:)'); shading flat
colorbar
set(gca,'ylim',[-5 -2],'xtick',0:6:36,'xticklabel',[0:6:18 0:6:12],'tickdir','out')
subplot(2,1,2)
pcolor(0:36,ledges,hdlep([1:end 1:end/2+1],:)'); shading flat
colorbar
set(gca,'ylim',[-5 -2],'xtick',0:6:36,'xticklabel',[0:6:18 0:6:12],'tickdir','out')

% linear scale diurnal histogram
subplot(2,1,1)
pcolor(0:36,edges,log2(hep([1:end 1:end/2+1],:))'); shading flat
colorbar
hold on
plot(-0.5:36.5,epsm([end 1:end 1:end/2+1],2),'g')
set(gca,'ylim',[0 2e-3],'xtick',0:6:36,'xticklabel',[0:6:18 0:6:12],'tickdir','out')
subplot(2,1,2)
pcolor(0:36,edges-2e-3,hdep([1:end 1:end/2+1],:)'); shading flat
colorbar
set(gca,'ylim',[-0.5e-3 2e-3],'xtick',0:6:36,'xticklabel',[0:6:18 0:6:12],'tickdir','out')

%% log2-frequency diurnal histogram
clf
pcolor(0:36,edges*1e3,log2(hep([1:end 1:end/2+1],:))'); shading flat
a=b2rcolormap(20);
caxis([-0.5 5])
cm=colormap([[1 1 1]; a([11 11:19],:)]);
hc=colorbar;
set(hc,'ylim',[-0.1 4.6],'ytick',[0 1:0.5:4.5],'yticklabel',[1 2 3 4 6 8 12 16 23] ,'fontsize',14)
hold on
plot(mod(yd*24-5,24),ep*1e3,'k.','markersize',3)
plot(24+mod(yd*24-5,24),ep*1e3,'k.','markersize',3)
hml=plot(-0.5:36.5,epsm([end 1:end 1:end/2+1],:)*1e3,'color',[0 0 0]);
set(hml(2),'linewidth',2)
plot(0.5:35.5,epse([1:end 1:end/2])*1e3,'ko') % local hour
set(gca,'ylim',[0 2],'xlim',[0 36],'xtick',0:6:36,'xticklabel',[0:6:18 0:6:12],...
    'tickdir','out','fontsize',16)
xlabel('local hour')



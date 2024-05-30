function [sthr,meannoise,kthr]=HildebrandSekhon3(S,z1side)
% [sthr,meannoise,kthr]=HildebrandSehkon(F,S)
% White noise threshold and mean finder modified from (Hildebrand and Sekhon 1974) 
% for finding noise beside inertial spectra. Thresholding is along
% frequency dimension, rather than spectral power. R2_cond=2 recommended.
% 2018-05-30 Simon de Szoeke

% Basic idea: exclude large spectral estimates until the remainder
% cannot be distinguished from white noise.
% My method includes lower and lower frequency spectral estimates until one can be
% distinguished from white noise.

% 2018-07-09 SPdeS: found a bug in calculation of variance Q. Fixing it
% makes the noise threshold too conservative, i.e. excludes too much of the
% signal. An alternative is to find the end of white noise by assuming it
% is normal and testing for when spectral estimates below a wavenumber
% exceed a value that makes them very unlikely to belong to the
% noise distribution.

%p1side=0.98;
%z1side=norminv(p1side);

ns=length(S)-1; % ignore Nyquist spectral estimate
%[s,ord]=sort(S(1:ns)); % sort s in order of power
ord=(ns:-1:1)';  % order from high frequency to low frequency
s=S(ord);
reord=zeros(size(s));
reord(ord)=1:ns; % indices reverse ordering so s(reord) == S(1:ns)

% Divide each sum by the length of the sum.
% but this results in too conservative a condition, throwing out too many
P=cumsum(s)./(1:ns)'; % <s> for increasingly larger ensembles of s for lower frequencies
P2=P.*P;
% Q=<s'^2> for increasingly larger ensembles of larger/low frequency s
Q=cumsum(s.*s)./(1:ns)' - P2; % correct HS eqn 7, SPdeS 2018-07-09
%{
R_2=P2./Q;
% might have to retune R_2cond for this correction; it appears to be even
% more conservative once corrected. R_2cond ~>=10 might be good.

% oops, this line isn't equivalent to HS eqn 7: SPdeS 2018-07-09
%Q=cumsum((s-P).^2)./(1:ns)'; % <s'^2> for increasingly larger ensembles of larger/low frequency s
%R_2=P2./Q;

% Standard noise threshold is when R_2 decreases to 1.
% Collections of spectral estimates S with R_2>1 are white noise;
% collections with R_2<1 are probably signal.
% s(1:ithr-1) are small noise
ithr=20-1+find(R_2(20:end)<R_2cond,1,'first'); % index of highest-wavenumber nonnoise spectral estimate
%}

% Better threshold:
% Find index of SE left (at lower wavenumber) of which 
% all standardized SEs exceed z1side. i.e. where the SEs first drop into
% the noise distribution.
firstnoise=find((s(reord)-P(reord))./sqrt(Q(reord))<z1side,1,'first'); % indexes original order S
% P, Q are mean and var of the distribution, s is its lowest-wavenumber realization
kthr=firstnoise-1; % last nonnoise SE index, S(kthr)>noise but S(kthr+1)<noise.
if isempty(kthr) || kthr==0
    kthr=1; % the spectrum is all white noise
end
sthr=S(kthr);
% it would be conservative to exclude kthr from inertial spectrum.

ithr=ord(kthr); % s(ithr) is the highest-wavenumber nonnoise spectral estimate
% average noise
meannoise=mednmean(s(1:ithr-4),5); % mean of 5 middle points, insensitive to outliers, smoother than median

%{
% Find the probability that the first nonnoise spectral estimate belongs to
% the empirically determined distribution of white noise spectral estimates.
ncd=normcdf(s,P,sqrt(Q));
% lastnoise=ns-find(ncd(ord)<p1side,1,'first'); % search backwards to avoid false positives

% try alternate threshold, modeling white noise as a small sample from a
% a normal distribution
t=(s([2:ns ns])-P)./sqrt(Q); % statistic for the next (lower frequency) spectral estimate
% posterior=p1side.^(1:ns);
% ttest=tinv(posterior,1:ns);
ttest=tinv(p1side,(1:ns)');
t2=(s([2:ns ns])-P)./sqrt(Q./(1:ns)'); % statistic for the next (lower frequency) spectral estimate
ttest2=tinv(0.98,1); % play with degrees of freedom: What do you do when trying to compare a single realization to a finite sample?
% jthr=find(t>ttest,1,'first');
jthr=min(ns,find(t<ttest,1,'last')+2); % more robust condition for emerging permanently above noise
% +1 bc we tested for the last that failed, and +1 bc t is advanced 1

% % Bayesian probability estimate --doesn't work
% beta=0.1; % adjustable design parameter: ratio of P(noise|S>r)/P(~noise|S>r)
% quotient=(1:ns)'./(ns-(1:ns)');
% pqrratio=tcdf(t,(1:ns)')*0.5.*quotient;
% % pqrratio, likelihood of noise increases as n increases, because quotient
% % increases. Suspect quotient is wrong estimate of P(noise)/P(not noise).
% % But we expect spectral estimates to rise above noise for lower
% % wavenumbers.
%}

% Note: HB1974's white noise variance of *frequency* is sigma_N^2 = F^2/12, where
% F=Ndf is the frequency span resolved by the spectrum. This variance
% follows from the spectral density-weighted integral of f'^2 for a flat
% distribution for [0 F]. The span F grows as we include more specral
% estimates. Here we use a nondimensional set of frequencies f.

%{
%test plots
clf
loglog(F,S,'.-')
hold on
plot(F(16),S(16),'x')
plot(F([1 end]),meannoise*[1 1])
plot(F([1 end]),sthr*[1 1])
plot(F([1 end]),sj*[1 1])
%}

% commented out to save time, this calculation helps justify the R_1, R_2
% threshold
%{
f=(1:ns)'; % in order
sigN2=f.*f/12; % array increasing in frequency span
css=cumsum(s);
sig2=cumsum(s.*f(ord).*f(ord))./css - (cumsum(s.*f(ord))./css).^2; % HS1974 eqn 4
R_1=sigN2./sig2; % eqn 8

% s2=s.*s;
%}

% R_2cond=1 is too conservative for discerning white
% noise from k^-5/3 turbulence. It cuts off too much turbulence, and sets
% the mean noise too high.

%{
% like HS1974 Fig 2b, but with x-axis reversed
clf
semilogx(s/sthr,R_1,'.-')
semilogx(s/sthr,R_2,'.-')
hold on
plot(1,R_2(ithr),'x')
plot(s([1 end])/sthr,[1 1])
text(0.5,1,'noise')
text(2,1,'signal -->')

factr=1.3889; F53=F.^(5/3);
vls=factr*(2*pi/Utop)^(2/3)*(F53.*(S-meannoise)); % dissipation ^2/3
epsilon=vls.^1.5; % dissipation
plot(F,epsilon)
plot(F([4 16]),epsilon([4 16]),'o')
plot(F(5:kthr),epsilon(5:kthr),'linewidth',2)
%}

% It looks like the spectrum rolls off at f>F(16)=0.4375 Hz

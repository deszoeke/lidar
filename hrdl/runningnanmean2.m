function [lowpass, nobs] = runningnanmean2(timeseries,window,minnobs,dim,iter)
% [lowpass, nobs] = runningnanmean(timeseries,window,minnobs,dim,iter)
% returns the running mean of timeseries, with a symmetric
% boxcar window of length window. Running means with less than
% minnobs observations are returned as NaNs. Operates running mean
% along dimension dim; iterates running mean iter times.
% 
% (c) 2014 Dec 17 Simon de Szoeke 
    
% set defaults
if ~exist('iter','var')
    iter=1;
    if ~exist('dim','var')
        dim=1;
    end
end
      
halfwidth = floor((window-1)/2);
if minnobs<1
    minnobs=2*halfwidth*minnobs;
end
    
% more efficient to use an accumulator, and handles more dimensions
% operate running mean along dimension dim
sz0=size(timeseries);
ndim=length(sz0);
otherdims=1:ndim;
otherdims=otherdims(otherdims~=dim);
nt=sz0(dim);
nother=prod(sz0(otherdims));

% pad x with 1 column each at the begnning and end
x=[NaN(nother,1) reshape(permute(timeseries,[otherdims dim]),[nother nt]) NaN(nother,1)];
% x(neverythinelse,1+ntimes+1): operate along 2nd dimension of input matrix x

for repeat=1:iter
    isf=isfinite(x);
    isf(:,[1 end])=false;
    ofs=x(find(isf,1));
    if isempty(ofs)
        warning('Data is all NaN on iter=%i', repeat);
        y = x;
        nwin = zeros(size(x));
        break
    end
    x=x-ofs; % improve precision of sums
    x(~isf)=0;

    % uses difference of cumsums to compute running mean
    c=cumsum(x,2);
    n=cumsum(isf,2);
    sh1=min(1+(1:nt)+halfwidth  , nt+2); % shift from forward, let end rest on zero padding
    sh2=max(1+(1:nt)-halfwidth-1, 1   );
    nwin=n(:,sh1)-n(:,sh2);
    y=(c(:,sh1)-c(:,sh2))./nwin+ofs;
    y(nwin<minnobs)=NaN;
    x(:,2:nt+1)=y; % preps for next iteration
%     x(:,[1 nt+2])=NaN;
end

lowpass=ipermute(reshape(y,[sz0(otherdims) nt]),[otherdims dim]);
nobs=ipermute(reshape(nwin,[sz0(otherdims) nt]),[otherdims dim]); % final iteration only

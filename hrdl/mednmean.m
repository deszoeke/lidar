function m = mednmean(s,n)
% mednmean(s,n) returns the n-point mean of the middle n points in (sorted) s.
% mednmean(s,1) is the median
% mednmean(s,length(s)) is the mean, but runs slower than mean.
% 
% (c) Simon de Szoeke 2018-06-06

sx=sort(s(:));
nx=find(isfinite(sx),1,'last');
i1=[(nx+1)/2-(n-1)/2, (nx+1)/2+(n-1)/2];
ii=[max(1,floor(i1(1))):min(nx,floor(i1(2))) max(1,ceil(i1(1))):min(nx,ceil(i1(2)))];
m=mean(sx(ii));


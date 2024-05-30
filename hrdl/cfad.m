function a = cfad( w, edges )
%     wcfad(height, bins) = cfad( w, edges )
% Compute counts of w(time,height) within bins defined by edges,
% with time as the sampling dimension.
%     edges(n) <= bins(n) < edges(n+1)
% e.g.
%     edges = [-Inf -7:0.1:7 Inf];

% if ~exist(edges, 'var')
%     edges = [-Inf -7:0.1:7 Inf];
% end

[~, nk] = size(w);
ne = length(edges);
a = zeros( nk, ne-1 ); % ne-1 bins if edges include -+Inf

for k = 1:nk
    isf = isfinite(w(:,k));
    di = discretize( w(isf,k), edges );
    for i = 1:length(di)
        a(k, di(i)) = a(k, di(i)) + 1;
    end
end
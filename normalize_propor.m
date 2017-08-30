function [dataNorm] = normalize_propor(data)

%find the largest max-min and compute scaling coefficient
[n,d]=size(data);
ma=zeros(1,d);
mi=zeros(1,2);

for dim=1:d
    inds = find(~isnan(data(:,dim)) & isfinite(data(:,dim)));
    ma(dim)=max(data(inds,dim));
    mi(dim)=min(data(inds,dim));
end
% max difference between max and min of particular variable
r = max(ma-mi);

dataNorm = (data - repmat(mi,n,1)) / r;
function [sems,mus,stds] = nansem( x,dim )
%SEMS=NANSEM(X,[DIM]); standard errors of the mean with missing data (NaN)
%   compute standard errors of the mean of x along dimension dim
% Jacek P. Dmochowski, 2019
% send bug reports (or props) to dmochowski AT gmail DOT com
%
% inputs
%   x: the data, can be arranged however you want
%   dim: the dimension over which to compute stats
% outputs
%   sems: standard errors of the mean at each dimension
%   mus: the means at each dimension
%   stds: standard deviations of each dimension
if nargin<2, dim=ndims(x); end
%lendim=size(x,dim);
lendim=sum(~isnan(x),dim);
sems=nanstd(x,[],dim)./sqrt(lendim);
stds=nanstd(x,[],dim);
mus=nanmean(x,dim);  % get means while at it
end


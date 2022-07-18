%% function se = standardError(ds, dim) 
% Compute standard error of the mean: SD/sqrt(N)
% 
% Inputs: 
% - ds: a matrix of data. This can be any size. 
% - dim: the dimension over which to take the SEM over. The default is the last 
%   dimension. For instance, the last dimension could be subjects in a within-subjects experiment. 
% 
% Outputs: 
% - SEM: the standard error of the mean: the standard deviation divided by
%   the square root of the number of measurements. NaNs are treated as
%   missing values and ignored. 
% 
% by Alex White, 2018, at the University of Washington. 
% 

function SEM = standardError(ds, dim, weights) 

if nargin<2 || ~exist('dim','var')
    dim = ndims(ds);
end
%if ds is a vector, make sure we take SEM over the one dimension that
%matters:
if isvector(ds)
    dim = find(size(ds)>1);
end

if nargin<3 || ~exist('weights','var')
    weights = ones(1,size(ds, ndims(ds)));
end
weights = weights/sum(weights);  

%N: count how many non-nan measurements there are 
N = sum(~isnan(ds),dim);
try
    SEM = sqrt(nanvar(ds,weights,dim))./sqrt(N);
catch
    keyboard
end
%% function SEM = standardError(ds, dim, weights) 
% Compute standard error of the mean: SD/sqrt(N)
% 
% Inputs: 
% - ds: a matrix of data. This can be any size. 
% - dim: the dimension over which to take the SEM over. The default is the last 
%   dimension. For instance, the last dimension could be subjects in a within-subjects experiment. 
%    If omitted, the default dimension is the last in ds. 
% - weights: optional vector of nonzero positive weights to apply in
% computation of variance 
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
    if length(ds)>1
        dim = find(size(ds)>1);
    else
        dim = 1;
    end
end

if nargin<3 || ~exist('weights','var')
    weights = ones(1,size(ds, ndims(ds)));
elseif ~isvector(weights)
    error('weights must be a vector');
end
    
weights = weights/sum(weights);  

%N: count how many non-nan measurements there are 
N = sum(~isnan(ds),dim);
try
    SEM = sqrt(var(ds,weights,dim,'omitnan'))./sqrt(N);
catch
    keyboard
end
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
    % Determine which dimension to use--the last one
    dim = find(size(ds)~=1, 1, 'last');
    if isempty(dim), dim = 1; end
end
%if ds is a vector, make sure we take SEM over the one dimension that
%matters:
if isvector(ds)
    dim = find(size(ds)>1);
end

if nargin<3 || ~exist('weights','var')
    weights = ones(1,size(ds, dim));
elseif ~isvector(weights)
    error('weights must be a vector');
end

if length(weights)~=size(ds, dim)
    error('The length of weights must equal the size of input ds over dimension dim');
else
    if isvector(weights) && isvector(ds) 
        if all(size(weights') == size(ds)) %if they are just vectors to be transposed
            weights = weights';
        end
    end
end
    
%set weights to 0 where the data are NaN
if isvector(ds)
    weights(isnan(ds)) = 0;
else %reduce each weight proportion to how much of ds is NaN
    %first replicate w so it is same size as ds
    reshapeDims = ones(1,ndims(ds));
    reshapeDims(dim) = length(weights);
    wr = reshape(weights, reshapeDims);
    repDims = size(ds);
    repDims(dim) = 1;
    w = repmat(wr, repDims);
    
    %now set weights to 0 where corresponding element of ds is NaN
    w(isnan(ds)) = 0;
    %now sum over other dimensions to get weights back to a vector
    dimsToSum = setdiff(1:ndims(ds), dim);
    for di=dimsToSum
        w = sum(w, di);
    end
    weights = squeeze(w);
end

weights = weights/sum(weights);  

%N: count how many non-nan measurements there are 
N = sum(~isnan(ds),dim);
try
    SEM = sqrt(var(ds,weights,dim,'omitnan'))./sqrt(N);
catch
    keyboard
end
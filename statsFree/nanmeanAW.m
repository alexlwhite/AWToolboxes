function avg = nanmeanAW(x,dim)
% Take the average, excluding NaNs
% inputs: 
% - x: data to average. Can be vector or matrix. If x is a matrix, the
% output is averaged over rows, by default.  
% - dim: dimension to average over (optional)
% 
% output: 
% - avg: the mean, ingoring nans

%default dimension is 1
if nargin == 1
    if isvector(x)
        dim = find(size(x)>1);
    else
        dim = 1; 
    end
end

% Find NaNs and set them to zero
nans = isnan(x);
x(nans) = 0;

% Count up non-NaNs.
n = sum(~nans,dim);
n(n==0) = NaN; % prevent error of dividing by zerio
% Sum up non-NaNs, and divide by the number of non-NaNs.
avg = sum(x,dim) ./ n;


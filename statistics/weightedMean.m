function y = weightedMean(x,w,dim)
%weightedMean   Weighted Average or mean value.
%   For vectors, weightedMean(X,W) is the weighted mean value of the elements in X
%   using non-negative weights W. For matrices, weightedMean(X,W) is a row vector 
%   containing the weighted mean value of each column.  For N-D arrays, 
%   weightedMean(X,W) is the weighted mean value of the elements along the first 
%   non-singleton dimension of X.
%
%   Each element of X requires a corresponding weight, and hence the size 
%   of W must match that of X.
%
%   weightedMean(X,W,DIM) takes the weighted mean along the dimension DIM of X. 
%
%   Class support for inputs X and W:
%      float: double, single
%
%   Example:
%       x = rand(5,2);
%       w = rand(5,2);
%       weightedMean(x,w)
%
% wmean written by Adam Auton in  2009. 
% Modified by Alex White, 2022

if nargin<2
    error('Not enough input arguments.');
end


if nargin==2 
  % Determine which dimension SUM will use
  dim = find(size(x)~=1, 1, 'last');
  if isempty(dim), dim = 1; end
end


% Check that dimensions of X match those of W.
if(~isequal(size(x), size(w)))
    if isvector(w) && isvector(x) 
        if all(size(w') == size(x)) %if they are just vectors to be transposed
            w = w';
        else
            error('Inputs x and w must be the same size.');
        end 
    elseif isvector(w) %replicate it to be the same size of x
        reshapeDims = ones(1,ndims(x));
        reshapeDims(dim) = length(w);
        wr = reshape(w, reshapeDims);
        repDims = size(x); 
        repDims(dim) = 1;
        w = repmat(wr, repDims);
        %fprintf(1,'\n(%s) Reshaping weight vector w to be a matrix of same size as x\n',mfilename);
    else
        keyboard
        error('Inputs x and w must be the same size.');
    end
end

%Deal with NaNs in data x that may not match 0s in w
if all(size(x)==size(w))
    w(isnan(x)) = 0; %set weights to 0 where x has missing values
else
    error('failure to set inputs x and w to be same size');
end

% Check that all of W are non-negative.
if (any(w(:)<0))
    error('All weights, W, must be non-negative.');
end

% Check that there is at least one non-zero weight.
if (all(w(:)==0))
    error('At least one weight must be non-zero.');
end


y = nansum(w.*x,dim)./nansum(w,dim);
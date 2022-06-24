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

% wmean written by Adam Auton in  2009. 
% Modified by Alex White, 2022

if nargin<2
    error('Not enough input arguments.');
end

% Check that dimensions of X match those of W.
if(~isequal(size(x), size(w)))
    if isvector(w) && all(size(w') == size(x)) %if they are just vectors to be transposed
        w = w';
    else
        keyboard
        error('Inputs x and w must be the same size.');
    end
end

% Check that all of W are non-negative.
if (any(w(:)<0))
    error('All weights, W, must be non-negative.');
end

% Check that there is at least one non-zero weight.
if (all(w(:)==0))
    error('At least one weight must be non-zero.');
end

if nargin==2 
  % Determine which dimension SUM will use
  dim = find(size(x)~=1, 1, 'last');
  if isempty(dim), dim = 1; end
end


y = nansum(w.*x,dim)./nansum(w,dim);
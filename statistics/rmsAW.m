%%  function r = rmsAW(x)
% Alex White's simple version of the RMS (root mean squared) function . 
% Inputs: 
% - x, a vector (or matrix) 
% - dim: optional, the dimension over which to compute the mean. If x is a
% vector, dim defaults to the long dimension of the vector. Otherwise it
% defaults to dim=1. 
% 
% Output: 
% - r, the rms of x 

function r = rmsAW(x, dim)

if nargin<2
    if isvector(x)
        sz = size(x);
        dim = find(sz>1);
    else
        dim = 1;
    end

end

r = sqrt(mean(x.*x,dim));
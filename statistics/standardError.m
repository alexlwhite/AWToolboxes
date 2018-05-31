%% function se = standardError(ds) 
% Compute standard error of the mean: SD/sqrt(N)
% 
% Inputs: 
% ds: a matrix of data. This can be any size, but the last dimension must
% be the variable you want to take the SEM over. For instance, the last
% dimension could be subjects in a within-subject psychology experiment. 
% 
% Outputs: 
% SEM: the standard error of the mean: the standard deviation divided by
% the square root of the number of measurements. 
% 
% by Alex White, 2018, at the University of Washington. 
% 

function SEM = standardError(ds) 

N = size(ds,ndims(ds));
SEM = std(ds,0,ndims(ds))/sqrt(N);
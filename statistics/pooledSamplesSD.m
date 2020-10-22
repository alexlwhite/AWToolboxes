%% function [pooledSD, pooledN, grandMean] = pooledSamplesSD(ms, ns, SDs)
% By Alex White, 2020, at Stanford University 
%
% Given the means, SDs, and number of measurements for K sets of data 
% (e.g., the across-trial means for K conditions in an experiment),
% compute what the SD (std dev.) would have been if you just concatenated
% or pooled all the data points across the K sets together, from the beginning. 
% This is useful for computing a t-statistic for the comparison of a set of
% conditions against another set of conditions. 
% Inputs
% - ms: a [j x K] matrix of means for K individual sets. If j>1, then
%   we assume there are j separate experiments you want to calculate this for. 
% - ns: a [j x K] matrix of the number of measurements (ie trials) in each set. 
% - SDs: a [j x K] matrix of Standard Deviations for each set. 
% 
% Output: the pooled SD, pooled N, and grand means. Each of those will be
% of size [j x 1]
%
function [pooledSD, pooledN, grandMean] = pooledSamplesSD(ms, ns, SDs)

%global mean
grandMean = sum(ms.*ns, 2)./sum(ns,2);

%pooled number of measurements across conditions 
pooledN = sum(ns,2); 

%pooled SD, across the K datasets or conditions 
pooledSD = sqrt( sum( (ns-1).*SDs.^2 + ns.*(ms-grandMean).^2 , 2) ./ (pooledN-1));

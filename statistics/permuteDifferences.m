% function [CI,pval,permDist] = permuteDifferences(myStatistic,x,nReps,CIrange)
%
% Permutation test on a sample of pairwise differences in the vector x. 
% Simulates repeating the experiment nReps times, assuming there is no true
% difference between the original two conditions. So on each repetition
% there's a 50% chance that each value in x flips sign. The null
% distribution of differences therefore has a mean 0. 
% 
% This function computes a confidence interval (CI) on the null distribution 
% statistic (function handle 'myStatistic'). 
%
% It also computes a p-value, using the function getBootPs. The two-tailed
% p-value is basically twice the proportion of the bootstrapped
% distribution greater than the value of myStatistic(x) on the original
% sample. 
%
% Inputs:
%    myStatistic:       a handle to a function that takes in a vector and
%                       returns a scalar (like 'mean' or 'std')
%
%   x                   sample vector for generating the statistic and its
%                       confidence intervals
%
%   nReps               number of sample-with-replacement iterations
%                       (default 2000)
%  
%   CIrange             confidence interval range (default 68.27)

%
% Outputs:     
%   CI:                 confidence interval
%   pval:               p-value: computed from  proportion of the null distribution
%                       more extreme than the measured value. But because
%                       the test is two-tailed, the p-value is actually
%                       twice that proportion. 
%   permDist:           vector of bootstrapped values of statistic   (1xnReps)
% 
% by Alex White at the University of Washington, 12/2018

function [CI,pval,permDist] = permuteDifferences(myStatistic,x,nReps,CIrange)
% Deal with defaults

if ~exist('BCFlag','var')
    BCFlag = 1;
end

if ~exist('CIrange','var')
    CIrange = 68.27;  %corresponds to +/- 1 s.d. for a normal distribution
end

if ~exist('nReps','var')
    nReps=2000;
end


permDist = NaN(nReps,1);
for pi=1:nReps
    flipsign = rand(size(x))<0.5; 
    xf = x;
    xf(flipsign) = -1*x(flipsign); 
    permDist(pi)=myStatistic(xf);
end

pval = getBootPs(permDist, myStatistic(x));


%% Compute CI 

Lo = (100-CIrange)/2;
Hi = (100+CIrange)/2;

CI(1) = prctile(permDist,Lo);
CI(2) = prctile(permDist,Hi);

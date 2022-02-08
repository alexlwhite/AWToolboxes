%% script for computing useful statistics on some vector of data points 
% add the BayesFactor toolbox to your path...
% get it here: https://klabhub.github.io/bayesFactor/
% then fill in the correct in this command: 
%addpath(genpath('/Users/alwhite/Dropbox/MatlabToolboxes/bayesFactor'));

%% generate fake data 
N = 15; %how many subjects
mu = 1; %true distribution mean
sd = 1.5; %true distribution SD

%create vector d of random Gaussian data with mu and sd:
d = randn(1, N)*sd + mu; 

%% open a txt file to save results to 
sf = fopen('sampleStats.txt', 'w');

%% report the sample mean and SD 
sampleMean = mean(d);
sampleSD = std(d);

fprintf(sf, 'Mean of %i samples = %.3f, SD = %.3f\n', N, sampleMean, sampleSD);

%% get a bootstrapped 95% confidence interval (CI) of the mean: 
%function to compute the "statistic" we want to get a CI of: 
statFunct = @mean;

%how many "repetitions" of the fake experiments: 
nReps = 1000;

%Confidence interval range:
CIrange = 95; 

%whether to "bias correct" the confidence interval
BCFlag = true; 

CI = boyntonBootstrap(statFunct, d, nReps,CIrange,BCFlag); 
fprintf(sf, '95%% bootstrapped CI = [%.3f %.3f]\n', CI(1), CI(2));

%if the CI excludes 0, we can reject the null hypothesis that the true mean
%is 0. 

%% t-test
[~, pval, ~, tst] = ttest(d);

fprintf(sf, 'T-test: t(%i) = %.3f, p = %.4f\n', tst.df, tst.tstat, pval);

%% Bayes factor: 
bayesF = bf.ttest(d);
fprintf(sf, 'Bayes Factor = %.3f\n', bayesF);

%if the BF is >3, we have some confidence that the "alternate hypothesis"
%(mean >0) %is more likely than the "null hypothesis"  (mean = 0). 
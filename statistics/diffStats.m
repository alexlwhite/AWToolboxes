%% function [tStat, bayesFactor, CI, sigStars] = diffStats(diffs, statsF, BCFlag, CIRange, nReps, compVal)
% run statistics on a vector of differences (or really any numbers). Common
% usage is to determine statistics on the mean difference between two
% experimental conditions. This function computes and prints the result of
% a t-test, Bayes Factor, and bootstrapped confidence interval. 
% 
% Inputs: 
% - diffs: a 1xN vector of measurments, typically within-subjects differences
%   between two conditions. N is the number of measurements. 
% - statsF: handle to a file where the results will be printed. If
%   statsF==1, results get printed to the command window. 
% - BCFlag: this is a Boolean, input to boyntonBootstrap, determines
%   whether the bootstrapped confidence interval is "bias corrected." Deaults
%   to true. 
% - CIRange: a single number that specifies therange of the confidence
% interval. Defaults to 95. 
% - nReps: number of bootstrap repetitions. Defults to 1000. 
% - compVal: value against which the mean of diffs is compared. Defaults
%   to 0. 
% 
% Outputs: 
% - tStat: a structure with outputs from Matlab's ttest function, with three variables attached: 
%     tstat, the t-value; pval, the p-value, and df, the degrees of freedom. 
% - bayesFactor: 1 number, the BF for a t-test, computed with the BayesFactor toolbox
%   (https://klabhub.github.io/bayesFactor/) 
% - CI: the bootstrapped confidence interval, computed with
%   boyntonBootstrap
% - sigStars: a character string of asterisks, if the difference is
%   "significant." Used for adding stars to a plot to indicate significance.
%   If the CI includes 0, then the effect is deemed not significant and
%   sigStars is empty ''. Otherwise, if tStat.pval<0.001, sigStars = '***',
%   if  tStat.pval<0.01, sigStars = '**', if tStat.pval<0.05, sigStars = '*'.


function [tStat, bayesFactor, CI, sigStars] = diffStats(diffs, statsF, BCFlag, CIRange, nReps, compVal)


if nargin<3
    BCFlag = true;
end
if nargin<4
    CIRange = 95;
end
if nargin<5
    nReps = 1000;
end
if nargin<6
    compVal = 0;
end
%compare to some value:     
diffs = diffs-compVal;

%% compute stats:
[~, tp, ~, tStat] = ttest(diffs);
tStat.pval = tp;
bayesFactor = bf.ttest(diffs);
CI = boyntonBootstrap(@mean, diffs, nReps,CIRange,BCFlag); 

%% print
if statsF>0
    % report the sample mean and SD
    sampleMean = mean(diffs);
    sampleSD = std(diffs);
    sampleSEM = standardError(diffs);
    
    fprintf(statsF, 'Mean of %i samples difference from %.1f = %.3f, SD = %.3f, SEM = %.3f\n', length(diffs), compVal, sampleMean, sampleSD, sampleSEM);
    
    %report confidence interval
    fprintf(statsF, '%i%% bootstrapped CI = [%.3f %.3f]\n', CIRange, CI(1), CI(2));
    
    %report t-test
    if tp>0.0001 %print usual decimel number if bigger than 0.0001
        fprintf(statsF, 't(%i) = %.2f, p = %.4f; ', tStat.df, tStat.tstat, tp);
    else %print base-10 exponent otherwise (e.g., 2e-6
        fprintf(statsF, 't(%i) = %.2f, p = %.2u; ', tStat.df, tStat.tstat, tp);
    end
    %report Bayes Factor
    fprintf(statsF, 'BF = %.3f\n', bayesFactor);
end

%% how "significant" is the result, for the purpose of putting asterisks 
isSig = all(CI<0) || all(CI>0);
if isSig
    if tStat.pval<0.001
        sigStars = '***';
    elseif tStat.pval<0.01
        sigStars = '**';
    else
        sigStars = '*';
    end
else
    sigStars = '';
end
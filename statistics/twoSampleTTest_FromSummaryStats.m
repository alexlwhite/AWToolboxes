%%[mdiff, tval, pval] = twoSampleTTest_FromSummaryStats(means, SEs, nTrials, contrastWeights)
% by Alex White, 2020, at Stanford University 
% 
%This function computes an independent (unpaired) t-test, given the means and SEs from several
%independent "conditions", each with a certain number of trials, for a certain contrast between some of those conditions. 
% This was written for the purpose of computing "contrasts" on beta weights from a GLM of fMRI data. 
%
% Inputs: 
% - means: a VxC matrix of across-trial mean values for each of C
%   conditions, for each of V "cases" (for example, MRI voxels); 
% - SEs: a VxC matrix of across-trial standard errors; 
% - nTrials: a 1xC vector of the number of trials in each condition 
% - contrastWeights: a 1xC vector of weights for the contrast. The positive
% weights must sum to 1, as must the negative weights. The t-test is
% for the contrast between the "positive" and "negative" conditions. 
% 
% This is challenging when there is more than 1 positive or negative
% condition. To compute the SD of the aggregated positive (and negative)
% condition, I use the function pooledSamplesSD, which computes
% what the SD would have been had you just but all the trials into 1 bin. 
% 
% according to Wikipedia, 
% t = (mean1 - mean2)/(sp * sqrt(1/n1+1/n2))
% sp = sqrt( ((n1-1)*var(x1) + (n2-1)*var(x2)) / (n1 + n2 - 2))
% https://en.wikipedia.org/wiki/Student%27s_t-test#Equal_or_unequal_sample_sizes,_similar_variances
% matlab's ttest2 uses the same equations. 
% 

function [mdiff, tval, pval] = twoSampleTTest_FromSummaryStats(means, SEs, nTrials, contrastWeights)

nV = size(means,1);

%pull out standard errors for each condition
posConds = find(contrastWeights>0);
negConds = find(contrastWeights<0);

wPos = sum(contrastWeights(posConds));
wNeg = sum(contrastWeights(negConds));
if abs(wPos-1)>10^-10 || abs(wNeg+1)>10^-10
    keyboard
    error('contrastWeights not set properly. The positive weights should sum to 1, as should the negative weights');
end

%mean difference between positive and negative conditions:
%Vx1  = VxC  *  Cx1
mdiff = means*contrastWeights;

%pull out standard errors from both positive and negative conditions
pSEs = SEs(:,posConds);
nSEs = SEs(:,negConds);

%get number of trials for each cond and repmat across voxels
pNTs = repmat(nTrials(posConds),nV,1);
nNTs = repmat(nTrials(negConds),nV,1);

%convert from standard errors to SDs
pSDs = pSEs.*sqrt(pNTs);
nSDs = nSEs.*sqrt(nNTs);

%get pooled SDs and number of trials over multiple conditions
if length(posConds)>1
    [pSDs, pNTs] = pooledSamplesSD(means(:,posConds), pNTs, pSDs);
end
if length(negConds)>1
    [nSDs, nNTs] = pooledSamplesSD(means(:,negConds), nNTs, nSDs);
end

%degrees of freedom for a two-sample, unpaired t-test
df = pNTs + nNTs - 2;

%pooled standard deviation (sp)
sp = sqrt( ( (pNTs - 1).*(pSDs.^2) + (nNTs - 1).*(nSDs.^2) ) ./ df  );

%standard error of the difference, the t-value denominator
seDiff = sp.*sqrt(1./pNTs + 1./nNTs);

%t-value is the ratio of difference in  means to standard error of that difference:
tval = mdiff./seDiff;

%two-tailed p-value
pval = 2*tcdf(-abs(tval),df);

% compare this to much simpler calculation suggested by KNK:
%simple approximation of the poolsed standard error, the denominator of the t-statistic:
%the sqrt of the sum of squared SEs from the 'active' and 'control'
%condition. This was suggested by KNK
%But first average over sub-conditions in each one:
%I think this only works for comparing 1 condition with 1 condition,
%when the numbers of trials are equal

%     contrast_se2 = sqrt(mean(pSEs,2).^2 + mean(nSEs,2).^2);
%     contrast_t2 = mdiff./contrast_se2;
%     figure; hold on;
%     plot(contrast_t, contrast_t2, '.');
%     xlabel('my t'); ylabel('KNK t');
%     axlims = [min([contrast_t; contrast_t2]) max([contrast_t; contrast_t2])];
%     plot(axlims, axlims, 'k-');
%     plot(axlims, 2*axlims, 'r-');
%     xlim(axlims); ylim(axlims); axis square;


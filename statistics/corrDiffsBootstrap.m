%% function [bootCI, bootP, meanCorrDiff] = corrDiffsBootstrap(pair1, pair2)
% 
% Bootstrap the difference between two correlations 
% Inputs
% - pair1: a Nx2 matrix of values. The idea is to compare the values in the
%  1st column against the second column, and compare that to how well the
%  columns in pair2 correlate with each other. 
%- pair2: another Nx2 matrix of values. 
% 
% outputs: 
% - bootCI: 95% bootstrapped confidence interval on the difference between
% the correlation coefficients between pair1 and between pair2
% - bootP: a simple p-value (the proportion of bootstrapped differences on
% the wrong side of 0). 
% - meanCorrDiff: mean of bootstrapped differences. 

function [bootCI, bootP, meanCorrDiff] = corrDiffsBootstrap(pair1, pair2)


nBoot = 100000;
N = size(pair1,1);
if size(pair2,1)~=N
    error('Pairs dont have same sample size');
end

cdiffs = NaN(nBoot,1);

for bi=1:nBoot
    is = randsample(1:N, N, true);
    cdiffs(bi) = corr(pair1(is,1), pair1(is,2)) - corr(pair2(is,1),pair2(is,2));
end

bootCI = prctile(cdiffs,[2.5 97.5]);
meanCorrDiff = mean(cdiffs);

%figure; histogram(cdiffs, 100);

if meanCorrDiff<0
    bootP = mean(cdiffs > 0); % Calculate the p-value based on the bootstrap distribution
else
    bootP = mean(cdiffs < 0); % Calculate the p-value based on the bootstrap distribution
end
disp(['Bootstrap p-value: ', num2str(bootP)]);

%% Gemini's justification for using bootstrapping and resampling data points with replacement
%% which means multiple copies of the same data points go into correlations: 

% Here is why drawing duplicate participants does not distort your correlation:
% 	You sample whole participants (rows), preserving structure: On each iteration, you draw an entire participant's set of scores [A_i,B_i,C_i,D_i ] together. Because the paired measurements within each participant are never broken apart, the internal covariance structure among conditions remains completely intact.
% 	Duplication alone does not change a correlation: Correlation measures relative bivariate alignment, not sample size. If you took a dataset of 20 people and cloned every single person 5 times to make 100 people, Pearson's r would remain identical.
% 	Variation creates the sampling distribution: Because a bootstrap sample randomly duplicates some participants and omits others, each resampled dataset gets a slightly different mix of high and low scatter. Calculating Δr=r_AB-r_CD on each of these 5,000 resampled datasets builds an empirical distribution that models the true sampling variability of your effect.
% 	The sample acts as the population: Bootstrapping treats your N participants as a mini-universe representing the broader population. Drawing N rows with replacement simulates what would happen if you ran your experiment 5,000 separate times with new random draws of N observers from that population.
% The method would only fail if you resampled cells independently (e.g., drawing a value for Condition A from Participant 3 and Condition B from Participant 8), which would destroy the within-subject design. As long as you resample entire participant rows intact, bootstrapping yields an unbiased test for comparing dependent correlations.
% 

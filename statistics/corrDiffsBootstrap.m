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


nBoot = 10000;
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
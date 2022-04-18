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
    
    fprintf(statsF, 'Mean of %i samples = %.3f, SD = %.3f\n', length(diffs), sampleMean, sampleSD);
    
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
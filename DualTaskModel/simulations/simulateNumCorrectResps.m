%% simulate the proportion of dual-task trials with 0, 1 or 2 correct responses

pProcessBoth = 0; %never process both stimuli
pStim1First = 0.5; %no bias to one side or the other

singleTaskDprimes = [1.7846   1.6386];
nDs = length(singleTaskDprimes);

serialModelAccs = NaN(nDs, 3);
fixedCapAccs = NaN(nDs, 3);

nSimTrials = 100000;
for di = 1:nDs

    %% serial model
    attnd1Mean = singleTaskDprimes(di);
    attnd2Mean = attnd1Mean;
    [~, ~, ~, ~, ~, dualTask_propNumCorrctResp] = SerialDualTaskModel(attnd1Mean, attnd2Mean, pStim1First, pProcessBoth, nSimTrials);
    serialModelAccs(di, :) = dualTask_propNumCorrctResp;

    %% fixed-capacity model 
    %assume balanced left/right focus
    pSamplesOnTask1 = 0.5;
    %dual-task dprime:
    dualD = singleTaskDprimes(di)./sqrt(1./pSamplesOnTask1);
    %if we assume neutral criterion, we can get prop correct: 
    neutralHR = 1-normcdf(dualD/2, dualD, 1); %hit ratewhen c is midway between 0 and d'
    neutralCR = normcdf(dualD/2, 0, 1); %correct reject rate 
    propCorr = mean([neutralHR neutralCR]); %proportion correct is mean of HR and CR

    %now if we have 2 of those coin flips independently, what's probability
    %of both, or just 1? 
    %could simulate 2 independent coin flips, or just calculate:
%         c1 = rand(1,nSimTrials)<propCorr;
%         c2 = rand(1,nSimTrials)<propCorr;
%         nCorrect = c1+c2;
%         fixedCapAccs(di,:) = [mean(nCorrect==0) mean(nCorrect==1) mean(nCorrect==2)];
       
   %calculate
   pBoth = propCorr^2;
   pNeither = (1-propCorr)^2;
   pOne = 1-pBoth-pNeither;
   fixedCapAccs(di,:) = [pNeither pOne pBoth];
end
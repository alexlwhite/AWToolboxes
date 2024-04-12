%% function [serialPred, fixedCapPred] = predictPropNumCorrectResps(singleDs, dualCrits, serialBias, parallelBias)
% This function predicts data in our dual-task experiments, specifically the proportion of dual-task trials with 0 correct responses, 1 correct response, or 2 correct responses. 
% It predicts that assuming the all-or-none, one-stim-per-trial serial
% model, and the fixed-capacity parallel model. 
% It takes as inputs the single-task d' levels for both sides (aka tasks),
% the estimated lexical decision criteria measured from dual-task trials,
% the estimated bias to process side 1 (task 1) according to the serial
% model, and the esimated bias to process side 1 according to the parallel
% model. 
% 
% Inputs: 
% - singleDs: a 1x2 vector of d' levels esimated for side 1 and side 2 in
%   the single-task (focused attention) condition 
% - dualCrits: a 1x2 matrix of SDT criterion levels for side 1 and 2 in the
%   dual-task condition. These are just norminv(FAR), where FAR is the false
%   alarm rate. In other words, the distance of the criterion from 0. 
%   We may also just assume the criterion is neutral. If so switch
%   assumeNeutralCriterion to "true" in the 1st line of the function. 
% - serialBias: a single number, in the range 0-1, that is the estimated proportion of dual-task trials when side 1 is
%   processed "first" (or rather, if the model assumes only 1 side is
%   processeed per trial, this is the proportion of trials when side 1 is
%   processed and side 2 is not). 
% - parallelBias: the esimated proportion of fixed capacity resources that
% are devoted to side 1, if that parallel model is true. 
% 
% Outputs: 
% - serialPred: a 1x3 vector that summarized the proportion of dual-task
%   trials when both responses are wrong, 1 response is correct, and both
%   responses are correct (in that order). These values sum to 1. 
% - fixedCapPred: just like serialPred, but as predicted by the
%   fixed-capacity parallel model. 
% 
% by Alex L. White, Barnard College, 2023 
%
% Note: this is not actually very useful, methinks. It kinda just captures
% that the serial model predicts fewer correct responses overall. 

function [serialPred, fixedCapPred] = predictPropNumCorrectResps(singleDs, dualCrits, serialBias, parallelBias)

assumeNeutralCriterion = false;
 attnd1Mean = singleDs(1);
    attnd2Mean = singleDs(2);

    pProcessBoth = 0; %assume total serial model, just 1 stim per trial always
    nSimTrials = 100000;
    [~, ~, ~, ~, ~, dualTask_propNumCorrctResp] = SerialDualTaskModel(attnd1Mean, attnd2Mean, serialBias, pProcessBoth, nSimTrials);
    serialPred = dualTask_propNumCorrctResp;

    %% fixed-capacity model 
    pSamplesOnTask1 = parallelBias;

    %dual-task dprimes
    fixedD1 = singleDs(1)./sqrt(1./pSamplesOnTask1);
    fixedD2 = singleDs(2)./sqrt(1./(1-pSamplesOnTask1));

    fixedDs = [fixedD1 fixedD2];
    propCorr = NaN(1,2);
    for side=1:2

        %use actual criterion measured, or assume neutral
        if assumeNeutralCriterion
            crit = fixedDs(side)/2;
        else
            crit = dualCrits(side);
        end

        hr = 1-normcdf(crit, fixedDs(side), 1); %hit ratewhen c is midway between 0 and d'
        cr = normcdf(crit, 0, 1); %correct reject rate 

        propCorr(side) = mean([hr cr]); %proportion correct is mean of HR and CR

    end
       
   %calculate
   pBoth = prod(propCorr);
   pNeither = prod(1-propCorr);
   pOne = 1-pBoth-pNeither;
   fixedCapPred = [pNeither pOne pBoth];
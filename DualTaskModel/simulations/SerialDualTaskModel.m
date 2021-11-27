%% [Ag, AgValsByIndex, dualTaskAccCorr, corrValsByIndex, tradeoffAgs] = SerialDualTaskModel(attndStim1Mean, attndStim2Mean, pStim1First, pProcessBoth, nT)
% This function simulates performance in a "dual-task" paradigm according to
% a serial processing model.
% It simulates many trials and constructs ROC curves to estimate accuracy
% for judging 2 stimuli in 3 conditions: single-task judging stim 1;
% single-task juding stim 2; and dual task (judging both stims).
% It also computes correlations between dual-task
% response accuracies.
%
% Congruency effects? not yet
%
% by Alex L. White, summer 2019 at the University of Washington.
%
%%%% Inputs:
%
%
%- attndStim1Mean: mean of target-present distribution of "target evidence" for
%  stimulus #1 when fully attended. If you have real data, this can be
%  estimated from single-task accuracy level for stimulus 1.
%
%- attndStim2Mean: mean of target-present distribution for stimulus 2 when
%  fully attended.
%
%- pStim1First: amount of bias in processing towards stimulus number 1 in dual-task trials.
%
%   For the serial model, we assume that on some proportion of dual-task trials (1-pProcessBoth),
%   only 1 stimulus is attended and the other is not processed at all. dualTaskAttnBias is
%   the proportion of those dual-task trials in which stim 1 is attended and stim 2 is ignored.
%
%- pProcessBoth: the probability of processing both
%   stimuli fully in dual-task trials. Must be between 0 and 1. If 0, we
%   get the all-or-none serial model; if 1, we get unlimited capacity
%   model.
%
% - nT: number of trials to simulate
%
%%%% Outputs
%
% - Ag: A matrix of areas under the ROC curve. Dimensions are explained by
%       output AgValsByIndex.
%      dim1 = single-task 1; single-task 2; dual-task
%      dim2 = stim 1; stim 2;
%      dim3 = all; incongruent, congruent
%
% - AgValsByIndex: structure with 1 field for each dimension of Ag. Each
%   field is a cell array that labels each row/column along that dimension.
%
% - dualTaskAccCorr: correlation between accuracy of 2 responses on
%      dual-task trials, divided into sets of trials by which sides had a
%      target
%
% - corrValsByIndex: structure that indicates what is in each element of dualTaskAccCorr
%
% - tradeoffAgs: AROC values in the dual-task condition, separated by whether the other
% side's response was incorrect or correct
%


function [Ag, AgValsByIndex, dualTaskAccCorr, corrValsByIndex, tradeoffAgs] = SerialDualTaskModel(attndStim1Mean, attndStim2Mean, pStim1First, pProcessBoth, nT)

%% parameters

%Number of trials simulated for each attention condition
if ~exist('nT','var')
    nT = 100000;
end

%number of possible responses the observer chooses between on each trial.
%a small number like 4 (which we actually use in an experiment)
%understimates AUC, relative to the true d'. So for accurate simulations it
%may be better to have a higher number
nRespLevs = 25;

%mean of distribution of E when no target is present
targAbstMean = 0;

%SD of distribution of E when attending fully to the stimulus
sigmaAttn =  1;
%SD of distribution of E when ignoring
sigmaIgnr =  1;
%Note: we could also cause guessing for ignored stimuli but setting this
%sigma to be very high, rather than setting targPresMeanIgnr to 0.

%probability of swapping sides by accident (one way to generate congruency
%effects; currently not used).
pSwapSides = 0;
swapSides = rand(nT, 1)<pSwapSides;

%% Set the overall probability of attending to side1 only and side2 only in dual-task situation
%probability of processing  both stimuli in dual-task trials.
%Under all-or-none serial model this should be 0. But we may assume that in
%some proportion of trials pProcessBoth > 0.
if pProcessBoth<0 || pProcessBoth>1
    error('\npProcessBoth must be between 0 and 1\n');
end

%Probability of processing only to stim 1 or stim 2, during dual-task
%trials. Controlled by pStim1First, which is the tendency to
%attend to stimulus 1.
pProcess1Only = (1-pProcessBoth)*pStim1First;
pProcess2Only = (1-pProcessBoth)*(1-pStim1First);


%% Target-present distributions of sensory evidence E
% 3x2 Matrix: 3 attention conditions x 2 stimuli

% Means:

%what should the mean be for "ignored" stimuli? lets set it to a "default" distribution midway between
%the target-present and target-absent distribution, for each side
evidenceMeanIgnored = mean([attndStim1Mean targAbstMean;
    attndStim2Mean targAbstMean], 2);

%orrr it's the target absent distribution
%evidenceMeanIgnored = [targAbstMean targAbstMean];


presMs = [attndStim1Mean evidenceMeanIgnored(1); ...  %single-task side 1
    evidenceMeanIgnored(2) attndStim2Mean; ... %single-task side 2
    attndStim1Mean attndStim2Mean];       %dual-task

%Standard deviations:
presSDs = [sigmaAttn sigmaIgnr; ...
    sigmaIgnr sigmaAttn; ...
    sigmaAttn sigmaAttn];



%% Target-absent distributions of sensory evidence E
%Means:
%abstMs = ones(size(presMs))*targAbstMean;

%target-absent means are only set to the 'correct' value of targAbstMean
%when that side is attended. Otherwise they're set to evidenceMeanIgnored.
abstMs = [targAbstMean evidenceMeanIgnored(1); %single-task side 1
    evidenceMeanIgnored(2) targAbstMean;   %single-task side 2
    targAbstMean  targAbstMean];   %dual-task


%Standard deviations (assumed to equal SDs for target-present)
abstSDs = presSDs;


%% Set decision criteria to be spaced out equally  (assume no bias)
mms = [presMs(:) abstMs(:)];
sss = [presSDs(:) abstSDs(:)];

Cs = generateBalancedCriteria(mms,sss,nRespLevs);


%% labels
cueLabs = {'single-task 1','single-task 2','dual-task'};
nCueConds = numel(cueLabs);

sideLabs = {'side1','side2'};

congLabs = {'all','incongruent','congruent'};
nCongConds = numel(congLabs);


AgValsByIndex.cue = cueLabs;
AgValsByIndex.side = sideLabs;
AgValsByIndex.congruency = congLabs;


%% loop through each condition and compute Ag (area under ROC curve)
nS = size(presMs,2); %number of simultaneous stimuli

Ag = NaN(nCueConds,nS,nCongConds);

for condi = 1:nCueConds
    %% set up variables for each trial
    
    % Whether target is present or absent in each stimulus on each trial
    pres = rand(nT,nS)<0.5;
    
    %set up means and std devs for each stimulus on each trial
    allMs = NaN(nT,nS);
    allSs = NaN(nT,nS);
    for si=1:nS
        allMs(pres(:,si)==1,si) = presMs(condi,si);
        allSs(pres(:,si)==1,si) = presSDs(condi,si);
        
        allMs(pres(:,si)==0,si) = abstMs(condi,si);
        allSs(pres(:,si)==0,si) = abstSDs(condi,si);
    end
    
    
    %% in dual-task condition, set which stimuli *dont* get processed, by adjusting means and SDs
    %Find subsets of trials when either only stim1, or only stim2 are processed,
    %and adjust SDs and means as needed
    if strcmp(cueLabs{condi},'dual-task')
        sws = rand(1,nT); %random numbers
        stimIgnored = zeros(size(sws));
        stimIgnored(sws<pProcess1Only) = 2; %trials when attend-1 only, stim 2 gets ignored
        stimIgnored(sws>=pProcess1Only & sws<(pProcess1Only+pProcess2Only)) = 1; %trials when attend-2 only, stim 1 gets ingored
        %if stimIgnored==0, that means both stimuli get processed and we
        %leave allMs and allSs as they are originaly set up (with both being processed).
        
        %set means for "ignored" stimuli
        allMs(stimIgnored==1,1) = evidenceMeanIgnored(1);
        allMs(stimIgnored==2,2) = evidenceMeanIgnored(2);
        
        allSs(stimIgnored==1,1) = sigmaIgnr;
        allSs(stimIgnored==2,2) = sigmaIgnr;
        
    end
    
    
    %% simulate psychophysical responses
    Rs = simulateDualTaskTrials(allMs,allSs,Cs,swapSides);
    
    
    %% compute area under ROC curves
    for si=1:nS
        %also compute cognruency effects
        for congi=1:3
            if strcmp(congLabs{congi},'all')
                congTrials = true(size(pres,1),1);
            elseif strcmp(congLabs{congi},'incongruent')
                congTrials = pres(:,1)~=pres(:,2);
            elseif strcmp(congLabs{congi},'congruent')
                congTrials = pres(:,1)==pres(:,2);
            end
            
            [hr,fr] = computeROCRates(pres(congTrials,si),Rs(congTrials,si),1:nRespLevs);
            Ag(condi,si,congi) = computeAROC(hr,fr);
            
        end
    end
    
    if strcmp(cueLabs{condi},'dual-task')
        %% Tradeoffs: Compute Ag conditional on whether the other side was responded to correctly or incorrectly
        respPres = Rs>(nRespLevs/2);
        respCorrect = respPres == pres;
        
        otherCorrects = [0 1];
        tradeoffAgs = NaN(1,length(otherCorrects));
        
        thisSidePres = [pres(:,1); pres(:,2)];
        thisSideResp = [Rs(:,1); Rs(:,2)];
        otherSideCorr = [respCorrect(:,2); respCorrect(:,1)];
        
        for otherCorrI = 1:length(otherCorrects)
            otherCorr = otherCorrects(otherCorrI);
            otherCorrTrials = otherSideCorr==otherCorr;
            
            [hr,fr] = computeROCRates(thisSidePres(otherCorrTrials),thisSideResp(otherCorrTrials),1:nRespLevs);
            tradeoffAgs(otherCorrI) = computeAROC(hr,fr);
        end
        
        %% Correlations: Compute correlations in accuracy for dual-task responses, depending on target presence on the two sides
        
        condLabs = {'any','neither','side1Only','side2Only','both'};
        
        nConds = length(condLabs);
        dualTaskAccCorr = NaN(1,nConds);
        for targPres = 1:nConds
            switch condLabs{targPres}
                case 'any'
                    subTrials = true(nT,1);
                case 'neither'
                    subTrials = ~pres(:,1) & ~pres(:,2);
                case 'side1Only'
                    subTrials = pres(:,1) & ~pres(:,2);
                case 'side2Only'
                    subTrials = ~pres(:,1) & pres(:,2);
                case 'both'
                    subTrials = pres(:,1) & pres(:,2);
            end
            
            rhos = corr(respCorrect(subTrials, :));
            dualTaskAccCorr(targPres) = rhos(2);
        end
        
        corrValsByIndex.targetsPresent = condLabs;
        
    end
end



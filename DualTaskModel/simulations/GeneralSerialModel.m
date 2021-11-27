%% [Ag, AgValsByIndex, dualTaskAccCorr, tradeoffAgs, corrValsByIndex, tradeoffValsByIndex] = GeneralDualTaskModel(attndStim1Mean, attndStim2Mean, pStim1First, pProcessBoth, nT)
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
% - tradeoffAgs: AROC values in the dual-task condition, separated by
% side, whether the other side had a target or not, and whether the other
% side's response was correct or not 
% 
% - tradeoffValsByIndex: structure that guides the indices for tradeoffAgs

%
%%% Issues: 
% - what to put in single-task rows of presMs and presSDs, for uncued
% stimuli? 
% - with small number of response levels, the computed Ag will always be
% understimated slightly becuase the whole curve is not filled out
%  

function [Ag, AgValsByIndex, dualTaskAccCorr, corrValsByIndex, tradeoffAgs, tradeoffValsByIndex] = GeneralSerialModel(attndStim1Mean, attndStim2Mean, pStim1First, pProcessBoth, nT)

%% parameters

%Number of trials simulated for each attention condition 
if ~exist('nT','var')
    nT = 100000; 
end

%number of possible responses the observer chooses between on each trial. 
%a small number like 4 (which we actually use in an experiment)
%understimates AUC, relative to the true d'. So for accurate simulations it
%may be better to have a higher number 
nRespLevs = 50; 

%make plots? 
plotROC = false;
plotAOC = false;
plotBars = false;
plotAccCorr = false;

%mean of distribution of E when no target is present
targAbstMean = 0; 

%mean of target-present distribution for either stimulus WHEN IGNORED
%This only really applies to the serial model. It's what forces the observer to guess 
%when asked about the stimulus that was not 'attended' or processed on that trial
%ACTUALLY now we use evidenceMeanIgnored, set midway between target-present and target-absent distributions for each side.  
%targPresMeanIgnr = targAbstMean;

%SD of distribution of E when attending fully to the stimulus
sigmaAttn =  1;
%SD of distribution of E when ignoring
sigmaIgnr =  1;
%Note: we could also cause guessing for ignored stimuli but setting this
%sigma to be very high, rather than setting targPresMeanIgnr to 0.

%probability of swapping sides by accident
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

%Note: does it matter what we put in presMs(1,2) and presMs(2,1)?
%Those are for the stimuli that are not pre- or post-cued in single-task
%condition. Same for presSDs.
%It may matter for congruency effects! 

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
    

    %% Compute ROC curves
    colrs = [0.3 0.8 0.6; 0.7 0.2 0.4];
    
    if plotROC
        figure; hold on;
        plot([0 1],[0 1],'k-');
    end
    hs = NaN(1,nS);
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
            if plotROC && congi==1
                hs(si)=plot(fr,hr,'.-','Color',colrs(si,:));
            end
        end
    end
    
    if plotROC
        xlabel('false alarm rate');
        ylabel('hit rate');
        axis square;
        title(cueLabs{condi});
        legend(hs,{'left','right'},'Location','SouthEast');
    end

    %% calculate correlation between accuracy of two responses...
    if strcmp(cueLabs{condi},'dual-task')
        [dualTaskAccCorr, tradeoffAgs, corrValsByIndex, tradeoffValsByIndex] = computeAccuracyContingencies(Rs,pres,nRespLevs);
    end
    
    
end

%% plot AOC
%re-format data
%- as: 2x2 matrix of area-under-ROC-curve measures.
%    rows = single-task, dual-task
%    columns = left side, right side
%- es: a 2x2 matrix of error bars for each of the conditions in as

if plotAOC
    congI = 1;
    as = NaN(2,2);
    as(1,1) = Ag(1,1,congI);
    as(1,2) = Ag(2,2,congI);
    as(2,:) = squeeze(Ag(3,:,congI));
    es = zeros(size(as));
    
    plotOpt.markSz = 12;
    plotOpt.doLegend = true;
    plotOpt.doAxLabels = true;
    
    plotOpt.axLineWidth = 1;
    plotOpt.datLineWidth = 1;
    
    edgeColor = [46 63 15]/255;
    dualFillColor = [1 1 1];
    plotOpt.edgeColors = edgeColor([1 1],:);
    plotOpt.fillColors = [edgeColor; dualFillColor];
    
    figure; hold on;
    
    %Note: this function uses older functions DualTaskSerialModel and
    %solveDualTaskSerialModel to fit dual-task data point...curious to see if
    %they agree
    
    [pGuessWhenIgnored, pAttendTask1] = plotAOCWithPredictions(as,es,plotOpt)
end

%% plot bars of Ag
if plotBars
    plotCuexSidexCongr(Ag,AgValsByIndex);
end

%% plot correlations and Ag contingent on other side 
if plotAccCorr
    plotDualTaskAccCorr(dualTaskAccCorr,corrValsByIndex.targetsPresent)
    
    plotAg_OtherSideContingent(tradeoffAgs,tradeoffValsByIndex)
end




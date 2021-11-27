%% [Ag, AgValsByIndex, dualTaskAccCorr, contingentAgs, corrValsByIndex, contingentAgValsByIndex] = GeneralDualTaskModel(modelType, attndStim1Mean, attndStim2Mean, pFlip2Side2, dualTaskAttnBias, extraCapacity, nT)
% This function simulates performance in a "dual-task" paradigm according to 
% either a serial or parallel model.
% It simulates many trials and constructs ROC curves to estimate accuracy
% for judging 2 stimuli in 3 conditions: single-task judging stim 1;
% single-task juding stim 2; and dual task (judging both stims).
% It also computes congruency effects and correlations between dual-task
% response accuracies. 
% In this version (April 25 2017), *congruency effects* are implemented by
% assuming that on some proportion of trials, regardless of instructions, the observer attends only
% to one stimulus (#1 or #2) and responds according to whatever evidence was perceived
% on that side, regardless of which side is asked about. This is controlled
% by the input parameter pFlip2Side2. 
%
% by Alex White, Spring 2017
% 
%%%% Inputs: 
%
%- modelType: 1=serial, 2=parallel 
%
%- attndStim1Mean: mean of target-present distribution of "target evidence" for
%  stimulus #1 when fully attended. If you have real data, this can be
%  estimated from single-task accuracy level for stimulus 1.
%
%- attndStim2Mean: mean of target-present distribution for stimulus 2 when 
%  fully attended. 
%
%- pFlip2Side2: the probability of attending only to side 2 and responding to
% both sides based on what was perceived there. This is to generate biased
% congruency effects. Ranges from -1 (always flip to side 1) to 1 (always
% flip to side 2)
% 
%- dualTaskAttnBias: amount of bias in processing towards stimulus number 1 in dual-task trials. 
%
%   For the serial model, we assume that on some proportion of dual-task trials (1-pDualBoth), 
%   only 1 stimulus is attended and the other is not processed at all. dualTaskAttnBias is
%   the proportion of those dual-task trials in which stim 1 is attended and stim 2 is ignored. 
%
%   For the parallel model, we assume that both stimuli are always
%   processed to some extent. In this case, dualTaskAttnBias is the
%   proportion of "sensory samples" allocated to stimulus 1. 
%
%- extraCapacity: relaxes the contraints of the strictest models. 
% 
%   In the serial model, this is probability of attending to both
%   stimuli equally in dual-task trials. Must be between 0 and 1. If 0, we
%   get the all-or-none serial model; if 1, we get unlimited capacity
%   model. 
%   
%   In the pararell model, this parameter changes how many total
%   "samples" are available for the two stimuli in dual-task condition. 
%    If extraCapacity=0, then we have the fixed-capacity modle: both stimuli 
%    share the same total numer of samples avaialble to just one in the 
%    single-task condition. 
%    If extraCapacity=1, then both stimuli get the full amount of samples
%    they each get in single-task condition, and we have unlimited capacity
%    parallel processing.
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
% - contingentAgs: AROC values in the dual-task condition, separated by
% side, whether the other side had a target or not, and whether the other
% side's response was correct or not 
% 
% - contingentAgValsByIndex: structure that guides the indices for contingentAgs

%
%%% Issues: 
% - what to put in single-task rows of presMs and presSDs, for uncued
% stimuli? 
% - with small number of response levels, the computed Ag will always be
% understimated slightly becuase the whole curve is not filled out
% 
% July 23, 2019: pFlip2Side2 is just one complicated case of an attention
% error that could cause congruency effects. A simpler case is that you
% just swap both sides, on some proportion of trials, but that would
% probably need to be different for single-task and dual-task conditions. 
% 
% Another source of congruency effects in the dual-task condition could be
% a bias to report the same answer for both sides (based on whichever one
% has more evidence? 

function [Ag, AgValsByIndex, dualTaskAccCorr, corrValsByIndex, contingentAgs, contingentAgValsByIndex] = GeneralDualTaskModel(modelType, attndStim1Mean, attndStim2Mean,pFlip2Side2, dualTaskAttnBias,extraCapacity, nT)

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


%% For serial model, set the overall probability of attending to side1 only and side2 only in dual-task situation
if modelType==1
    %probability of processing  both stimuli in dual-task trials. 
    %Under all-or-none serial model this should be 0. But we may assume that in
    %some proportion of trials pDualBoth > 0.
    pDualBoth = extraCapacity; 
    if pDualBoth<0 || pDualBoth>1
        error('\nFor serial model, extraCapacity (pDualBoth) must be between 0 and 1\n'); 
    end
    
    %Probability of attending only to stim 1 and stim 2, during dual-task
    %trials. Controlled by dualTaskAttnBias, which is the tendency to
    %attend to stimulus 1.
    if dualTaskAttnBias>=0
        pDualAttnd1 = (1-pDualBoth)*dualTaskAttnBias;
        pDualAttnd2 = (1-pDualBoth)*(1-dualTaskAttnBias);
    else
        pDualAttnd2 = (1-pDualBoth)*-1*dualTaskAttnBias;
        pDualAttnd1 = (1-pDualBoth)*(1+dualTaskAttnBias);
    end
end

%% Target-present distributions of sensory evidence E
% 3x2 Matrix: 3 attention conditions x 2 stimuli 

% Means: 

%what should the mean be for "ignored" stimuli? lets set it midway between
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
%% For parallel model, modulate SDs of distributions of evidence when dividing attention       
if modelType == 2  
    T = 1 + extraCapacity;
    if T<=0 || T>2
        error('For serial model, extraCapacity must be > -1 and <= 1\n'); 
    end
    pSampOnStim1 = dualTaskAttnBias.*T;
    pSampOnStim2 = (1-dualTaskAttnBias).*T;
    pSampOnStim1(pSampOnStim1>1)=1;
    pSampOnStim2(pSampOnStim2>1)=1;
    
    %don't allow the proportion of smaples exactly zero, or the sigma goes
    %to infinity and causes a crash 
    pSampOnStim1(pSampOnStim1==0) = 0.000001; 
    pSampOnStim2(pSampOnStim2==0) = 0.000001;


    sigma1Div = sqrt(1/pSampOnStim1)*sigmaAttn;
    sigma2Div = sqrt(1/pSampOnStim2)*sigmaAttn;

    presSDs(3,:) = [sigma1Div sigma2Div]; 
end

%% Target-absent distributions of sensory evidence E
%Means:
%abstMs = ones(size(presMs))*targAbstMean;

%target-absent means are only set to the 'correct' value of targAbstMean
%when that side is attended. Otherwise they're set to evidenceMeanIgnored. 
abstMs = [targAbstMean evidenceMeanIgnored(1); 
         evidenceMeanIgnored(2) targAbstMean; 
         targAbstMean  targAbstMean]; 
     
 
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
    

    %% If serial model, adjust means and SDs if in dual-task
    %Find subsets of trials when either only stim1, or only stim2 are attended,
    %and adjust SDs and means as needed 
    if strcmp(cueLabs{condi},'dual-task') && modelType==1
        sws = rand(1,nT); %random numbers
        switches = zeros(size(sws));
        switches(sws<pDualAttnd1)=1; %trials when attend-left only
        switches(sws>=pDualAttnd1 & sws<(pDualAttnd1+pDualAttnd2)) = 2; %trials when attend-right only
        
        %set means for "ignored" stimuli
        allMs(switches==2,1) = evidenceMeanIgnored(1);
        allMs(switches==1,2) = evidenceMeanIgnored(2);
        
        allSs(switches==2,1) = sigmaIgnr;
        allSs(switches==1,2) = sigmaIgnr;
        

    end
    
    %% implement accidental switching of attention  (for congruency effects)
    %if pFlip2Side2 is positive, then on that fraction of trials we process
    %only side 2 (flip2Side = 2). 
    if pFlip2Side2>=0  
        flip2Side = 2;
        pFlip = pFlip2Side2;
    %otherwise if it's negative, then on -1*pFlip2Side2 fraction of trials
    %we process only side 1 (flip2Side = 1);
    else
        flip2Side = 1; 
        pFlip = -1*pFlip2Side2;
    end
    
    %switchSides: says which side to switch attention to fully on each
    %trial:
    switchSides = zeros(nT,1);
    switchSides(rand(nT,1)<=pFlip) = flip2Side;
        
    %set the flip2side  to be attended fully on some switch trials, by re-setting 
    %means and sigmas of distributions of E. (it doesn't matter how we set
    %the other side's distributions, because it's never responded to).
    allMs(pres(:,flip2Side)==1 & switchSides>0,flip2Side) = presMs(flip2Side,flip2Side);
    allSs(pres(:,flip2Side)==1 & switchSides>0,flip2Side) = presSDs(flip2Side,flip2Side);
    
    allMs(pres(:,flip2Side)==0 & switchSides>0,flip2Side) = abstMs(flip2Side,flip2Side);
    allSs(pres(:,flip2Side)==0 & switchSides>0,flip2Side) = abstSDs(flip2Side,flip2Side);
    
        
    %% simulate psychophysical responses
    Rs = simulateDualTaskTrials(allMs,allSs,Cs,switchSides);
    
    %% calculate correlation between accuracy of two responses...
    if strcmp(cueLabs{condi},'dual-task')
        [dualTaskAccCorr, contingentAgs, corrValsByIndex, contingentAgValsByIndex] = computeAccuracyContingencies(Rs,pres,nRespLevs);
    end
    
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
    
    plotAg_OtherSideContingent(contingentAgs,contingentAgValsByIndex)
end




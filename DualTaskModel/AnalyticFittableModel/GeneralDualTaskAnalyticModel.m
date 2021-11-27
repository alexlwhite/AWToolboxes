%% [Ag, dualTaskAccCorr, AgValsByIndex] = GeneralDualTaskModel(modelType, attndStim1Mean, attndStim2Mean, pFlip2Side2, dualTaskAttnBias, extraCapacity)

    %July 7: SOMETHING IS TOTALLY BROKEN. NEED TO RESTRUCTURE. 
    % Tried to make this function compute ROCs not by simulating tons of
    % tirals but just by comptuing, given mu and sigma, probability of each
    % response. Should work. Conditions are kinda complicated though. 


% This function computes accuracy in a "dual-task" paradigm according to
% either a serial or parallel model.
% It computes the probability of each response given signal-detection theory parameters mean
% and SD of the sensory evidence distribution, and decision criteria.
% Then it computes ROC curves to estimate accuracy.
% It does so for each of  conditions: single-task judging stim 1;
% single-task juding stim 2; and dual task (judging both stims).
% It also computes congruency effects and correlations between dual-task
% response accuracies.
% In this version (July 7 2017), *congruency effects* are implemented by
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
% congruency effects.
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
%%%% Outputs
%
% - Ag: A matrix of areas under the ROC curve. Dimensions are explained by
%       output AgValsByIndex.
%      dim1 = single-task 1; single-task 2; dual-task
%      dim2 = stim 1; stim 2;
%      dim3 = all; incongruent, congruent
%
% - dualTaskAccCorr: correlation between accuracy of 2 responses on
%      dual-task trials
%
% - AgValsByIndex: structure with 1 field for each dimension of Ag. Each
%   field is a cell array that labels each row/column along that dimension.
%
%
%%% Issues:
% - what to put in single-task rows of presMs and presSDs, for uncued
% stimuli?
% - with small number of response levels, the computed Ag will always be
% understimated slightly becuase the whole curve is not filled out

function [Ag, AgValsByIndex, dualTaskAccCorr, contingentAgs, corrValsByIndex, contingentAgValsByIndex] = GeneralDualTaskAnalyticModel(modelType, attndStim1Mean, attndStim2Mean,pFlip2Side2, dualTaskAttnBias,extraCapacity)

%% parameters

%number of possible responses the observer chooses between on each trial
nRespLevs = 4;

%make plots?
plotROC = true;
plotAOC = true;
plotBars = true;
plotAccCorr = true;

%mean of distribution of E when no target is present
targAbstMean = 0;

%mean of evidence distribution for either stimulus WHEN IGNORED
%This only really applies to the serial model. It's what forces the observer to guess
%when asked about the stimulus that was not 'attended' or processed on that trial
ignoredMean = 0;

%SD of distribution of E when attending fully to the stimulus
attndSigma =  1;
%SD of distribution of E when ignoring
ignoredSigma =  1;
%Note: we could also cause guessing for ignored stimuli but setting this
%sigma to be very high, rather than setting targPresMeanIgnr to 0.


%% For serial model, set the overall probability of attending to left and right only in dual-task situation
if modelType==1
    %probability of attending to both stimuli in dual-task trials.
    %Under all-or-none serial model this is 0. But we may assume that in
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
presMs = [attndStim1Mean ignoredMean; ...  %single-task left
    ignoredMean attndStim2Mean; ... %single-task right
    attndStim1Mean attndStim2Mean];       %dual-task

%Standard deviations:
presSDs = [attndSigma ignoredSigma; ...
    ignoredSigma attndSigma; ...
    attndSigma attndSigma];

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
    
    sigma1Div = sqrt(1/pSampOnStim1)*attndSigma;
    sigma2Div = sqrt(1/pSampOnStim2)*attndSigma;
    
    presSDs(3,:) = [sigma1Div sigma2Div];
end

%% Target-absent distributions of sensory evidence E
%Means:
abstMs = ones(size(presMs))*targAbstMean;
%Standard deviations (assumed to equal SDs for target-present)
abstSDs = presSDs;


%% Set decision criteria to be spaced out equally  (assume no bias, across all trials of all conditions)
mms = [presMs(:) abstMs(:)];
sss = [presSDs(:) abstSDs(:)];

Cs = generateBalancedCriteria(mms,sss,nRespLevs);

%Exclude first criterion which is negative infinity, doesn't make sense
Cs = Cs(2:end);

%% labels
cueLabs = {'single-task 1','single-task 2','dual-task'};
nCueConds = numel(cueLabs);

sideLabs = {'left','right'};

congLabs = {'all','incongruent','congruent'};
nCongConds = numel(congLabs);


AgValsByIndex.cue = cueLabs;
AgValsByIndex.side = sideLabs;
AgValsByIndex.congruency = congLabs;


%% loop through each condition and compute Ag (area under ROC curve)
nS = size(presMs,2); %number of simultaneous stimuli

Ag = NaN(nCueConds,nS);


for condi = 1:nCueConds
    allMs = zeros(1,2);
    allSs = zeros(1,2);
    respProbs = NaN(2,2,2,2,nRespLevs);
    
    for stim1Pres = [0 1]
        for stim2Pres = [0 1]
            
            arePres = [stim1Pres stim2Pres];
            for stimI = [1 2]
                if arePres(stimI)
                    allMs(stimI) = presMs(condi,stimI);
                    allSs(stimI) = presSDs(condi,stimI);
                else
                    allMs(stimI) = abstMs(condi,stimI);
                    allSs(stimI) = abstSDs(condi,stimI);
                end
            end
            
            if strcmp(cueLabs{condi},'dual-task')
                stimsToProcess = [1 2];
                if modelType==1
                    attnStates=[1 2];
                else
                    attnStates = 1;
                end
            else
                stimsToProcess = condi;
                attnStates= 1;
            end
            
            for asi = attnStates
                if modelType==1 %implement 1-at-a-time serial
                    allMs(3-asi) = ignoredMean;
                    allSs(3-asi) = ignoredSigma;
                end
                
                for spi = stimsToProcess
                    cumPs = normcdf(Cs,allMs(spi),allSs(spi));
                    respProbs(stim1Pres+1,stim2Pres+1,asi,spi,:) = diff([0 cumPs]);

                end
            end
        end
    end
    
    %if serial switching model in dual-task situation, take weigthed sum
    %across the two attentional states
    if strcmp(cueLabs{condi},'dual-task') && modelType==2
        attnStateDim = 3;
        thesePs = squeeze(respProbs(condi,:,:,:,:,:));
        sz = ones(1,5);
        sz(attnStateDim)=2;
        weights = zeros(sz);
        weights(1) = dualTaskAttnBias;
        weights(2) = 1-dualTaskAttnBias;
        repSz = size(thesePs);
        repSz(attnStateDim) = 1;
        weightsR = repmat(weights,repSz);
        
        weightedPs = sum(thesePs.*weightsR,attnStateDim);
        
        respProbs(:,:,1,:,:) = weightedPs;
        
        %then just get ride of the attentional state dimiension
        respProbs = squeeze(respProbs(:,:,1,:,:));
    end
    
    %COMPUTE area under ROC for each stimulus in this condition
    for stimI = 1:2
        % - HR: a vector of length C, containing the hit rates compute at C
        % different criteria, from conservative to liberal. HR should therefore be
        % monotonically increasing. It's best if HR includes the most extreme
        % criteria that produce hit rates of 0 and 1, to anchor the curve at the
        % corners (even if that wasn't really possible for the subject to report in
        % the experiment). The corresponding function computeROCRates should do
        % that. But if the 1st element of HR is greater than 0, we'll concatenate it
        % with a 0, and if the last element of HR is less than 1, we'll concatenate
        % a 1 onto it.
        % - FR: a vector (also of length C) containing the corresponding false
        % alarm rates. Similarly, it should begin at 0  and to go 1. If not, those
        % value are concatenated.
        
        if stimI==1 %this is crazy
            presRespPs = squeeze(respProbs(2,:,1,:));
            abstRespPs = squeeze(respProbs(1,:,1,:));
        else
            presRespPs = squeeze(respProbs(:,2,2,:));
            abstRespPs = squeeze(respProbs(:,1,2,:));
        end
        %average over other side pres/abst (CRAZY)
        presRespPs = mean(presRespPs,1)
        abstRespPs = mean(abstRespPs,1)
        
        hr = cumsum(fliplr(squeeze(presRespPs)))
        fr = cumsum(fliplr(squeeze(abstRespPs)))
        Ag(condi,stimI) = computeAROC(hr,fr);
        
        
        if plotROC
            colrs = [0.3 0.8 0.6; 0.7 0.2 0.4];
            
            figure; hold on;
            plot([0 1],[0 1],'k-');
            hs(stimI)=plot(fr,hr,'.-','Color',colrs(stimI,:));
            
        end
        
        keyboard
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
    as(1,1) = Ag(1,1);
    as(1,2) = Ag(2,2);
    as(2,:) = squeeze(Ag(3,:));
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

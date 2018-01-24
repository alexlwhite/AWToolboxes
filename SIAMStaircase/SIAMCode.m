%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Bits of code to run the SIAM staircase %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%     Alex L White, December 2014        %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Set parameters, into a structure "stair"
stair.t                   = 0.5;   %desired maximum reduced hit rate (MRHR = HR-FAR). 0.5 is a good choice for intermediate performance. See Kaernbach (1990). 
stair.startStep           = 0.3;   %starting step size 
stair.nRevToStop          = 9;     %number of "good" reversals (with stable step size) to complete before terminating
stair.revsToHalfInt       = [1 2]; %on which reversals to halve step size, starting from first trial or starting just after step size reset
stair.revsToReset         = 100;   %on which reversals to reset, in case staircase continues
stair.nStuckToReset       = 5;     %the number of sequential hits all at the same intensity at which step size is reset
stair.threshType          = 1;     %How to estimate threshold: 1, for mean of reversal values; or 2, for all intensity values.

stair.terminateByReversalCount   = true;


%% Initialize staircase 
startLev = log10(0.3);          %starting level of staircase. If just 1 staircase, start at best guess of threshold. If 2 interleaved, start 1 high and 1 low. 
allLevs  = log10(0:0.001:1);    %all available stimulus levels 

stair.ss = initSIAM(stair.t, stair.startStep, startLev, allLevs, stair.revsToHalfInt, stair.revsToReset, stair.nStuckToReset); 

%% On each trial, extract recommended intensity level
stimLevel = stair.ss.intensity;

%% Then at the end of the trial, update the staircase based on what happened 
% targPres, whether the signal was present (1) or absent (0) on this last trial
% chosenRes, whether the subjects said yes (1) or no (0) 
stair.ss = updateSIAM(stair.ss,targPres,chosenRes);


%% Determine whether the staircase is completed 
if stair.terminateByReversalCount
    stairDone = stair.ss.nRevStableSteps>=stair.nRevToStop;
else
    stairTone = trialNum > totalTrials;
end

%% At the end of experiment, estimate threshold 
[thrsh, ntrls, sd] = estimateSIAM(stair.ss,stair.threshType);
%thrsh = threshold estimate 
%ntrials = number of trials with intensity values averaged to get thrsh
%sd: standard deviation of threshold estimate 

%% Plot staircase 
plotSIAM(stair.ss,stair.threshType);


%% Some other considerations: 
%- you might want to ignore a few trials at the start of a block before
%  using the staicase (and on those first trials just leave intensity at
%  starting level)
%- Often it's a good idea to control the stimulus in log units. That is, give the intensity to the staircase
% (initSiam and upDateSiam) in log units, and when you extract it, convert it back.
%function ss = initSIAM(t, startStep, startI, iSet, revsToHalveI, revToReset, nStuckToReset) 
%
% 2012 by Alex White
%
% initialize a SIAM staircase following Kaernbach, 1990
%
% Inputs:
% - t is the desired performance level, in units of the maximum reduced hit
% rate (hit rate - false alarm rate), which varies between 0 and 1. 0.5 is
% a good choice. 
% We assume that signal and noise trials are equally likely 
% - startStep is the starting intensity step size. 
% - startI is the starting intensity level 
% - iSet is the set of possible stimulus levels, if lenght(iSet)>2.
%   If iSet is a vector of length 2, it is taken instead as the minimum and 
%   maximum stimulus levels and the intensity level is otherwise free to wander. 
% - revsToHalveI is a vector of reversal numbers after which to halve the step size
% - revToReset is the number of reversals since last reset (or start) at which to reset the step size to the initial value
% - nStuckToReset is the number of sequential hits stuck at the same intensity level after which step size is reset to initial value 
%
% Outputs: 
% - ss, a siam staircase structure

function ss = initSIAM(t, startStep, startI, iSet, revsToHalveI, revToReset, nStuckToReset) 

ss.t             = t;
ss.stepSz        = startStep;
ss.intensity     = startI;

%iSet: if more than 2, its the set list of intensitiy levels to use. 
if length(iSet) > 2
    ss.iSet      = iSet;

%if only length 2, its the min and max of allowable levels
elseif length(iSet) == 2
    ss.iSet      = [];
else
%if empty, intensity levles are totally free     
    ss.iSet      = [];
end
ss.min           = min(iSet);
ss.max           = max(iSet);

ss.revsToHalveI  = revsToHalveI; 
ss.revToReset    = revToReset;
ss.nStuckToReset = nStuckToReset; 

%step size is computed for each type of response: hit, miss, false alarm,
%and correct rejection, in that order 
steps = [-1 t/(1-t) 1/(1-t) 0]; 
if any(steps(2:3)<1) 
    steps = steps/abs(min(steps(2:3)));
end
ss.steps = steps*ss.stepSz;
ss.origSteps = ss.steps;


ss.pres             = [];         %vector of 1s and 0s indicating whether signal was present on each trial
ss.resp             = [];         %vector of 1s and 0s indicating whether the subject said yes or no on each trial
ss.ints             = [];         %vector of stimulus intensities used on each trial
ss.reversalInts     = [];         %vector of intensity levels at reversals
ss.reversalTs       = [];         %vector of trial numbers of reversals
ss.revStableStepsTs = [];         %vector of trial numbers of reverals in which step size was not changed 
ss.revStableStepsIs = [];         %vector of intensitly levels at reverals in which step size was not changed 
ss.nreversals       = 0;          %count of number of reversals so far
ss.nRevStableSteps  = 0;          %count of number of reversals in which step size was not changing 
ss.resetT           = 1;          %trial numbers at which counter and step size reset to original
ss.resetRev         = 0;          %reversal number at which counter and step size reset
ss.tnum             = 0;          %count of trials done
ss.direction        = [];         %whether last step was up or down


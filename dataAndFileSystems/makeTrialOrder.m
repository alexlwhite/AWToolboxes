% trials = makeTrialOrder(design,numTrials)
% by Alex White
%
% This function assigns experimental parameters to a set of trials
% in a fully randomlized, counterbalanced fasion.
%
% Inputs:
% - design: a structure with fields "parameters" and "randVars."
%
%   parameters is a structure with fields corresponding to the experimental parameters
%   that you want to be fully counterbalanced. The name of each field is the name of the parameter. Each
%   field should be a vector containing the values that parameter can take.
%   For example, design.parameters.targetPosition=[1 2]; and design.parameters.targetOrientation=[-45 45]
%   would ensure that you have an equal number of trials with the target at
%   position 1 and orientation -45, at position 2 and orientation 45, etc.
%
%   randVars is a structure with fields corresponding to the
%   experimental parameters that you don't need fully counterbalanced. For
%   example, design.randVars.targetDelay=[0.2 0.25 0.3 0.35 0.45]
%   would be used to set a delay to one of those values, randomly on each
%   trial, but the function will not ensure that there are an equal number
%   of each on each trial.
%
% - numTrials: the total number of trials in
%   the block you desire. If the number of trials is not a multiple of the minimum
%   number for counterbalancing, this function will print a warning but
%   will return exactly how many trials you requested (even if that means
%   the counterbalancing is imperfect). 
%
% Output:
% - trials: a table, with one row for each trial, and one column for each variable 
%   (from both ''parameters'' and 'randVars'')
% - minCounterbalanceTs: minimum number of trials necessary for perfect counterbalancing

function [trials, minCounterbalanceTs] = makeTrialOrder(design, numTrials)

trials = table;
%%%%%%%%%%%%%%%%
%calculate how many trials to counterbalance everything
%%%%%%%%%%%%%%%%
if isfield(design, 'parameters')
    params=fieldnames(design.parameters);
    
    %compute how many trials it takes to fully counterbalance 
    nParams = length(params);
    paramLengths = NaN(1,nParams);
    for p=1:nParams
        eval(sprintf('thesePs = design.parameters.%s;', params{p}));
        paramLengths(p) = length(thesePs);
        %make them all into column vectors
        if size(thesePs,1)==1
            thesePs = thesePs';
            eval(sprintf('design.parameters.%s = thesePs;', params{p}));
        end
    end
    
    %minimum number of trials to get 1 trial of each possible combination
    %of all the parameters (i.e., for perfect counterbalancing)
    minCounterbalanceTs = prod(paramLengths);
    
    %number of repetitions of each trial type, given numTrials requested 
    numReps=ceil(numTrials/minCounterbalanceTs);
    
    %number of trials assuming all trial types occur equally often 
    fullNumTrials = numReps*minCounterbalanceTs;
    
    if numTrials<minCounterbalanceTs
        fprintf(1,'\n(makeTrialOrder) WARNING: Asked for doing %i trials, but %i required for 1 trial of each parameter combination.\n', numTrials, minCounterbalanceTs);
    end
    if mod(numTrials,minCounterbalanceTs)~=0
        fprintf(1,'\n(makeTrialOrder) WARNING: Asked for %i trials, which is not a multiple of the %i trials required for 1 of each parameter combination.\n', numTrials, minCounterbalanceTs);
    end
    
    %%%%%%%%%%%%%%
    % Compute Counterbalancing
    %%%%%%%%%%%%%%
    %% loop thru all possible combinations and add a trial
    bigCmd = ''; t=0;
    for p=1:nParams
        eval(sprintf('%s = NaN(minCounterbalanceTs,1);', params{p}));
        bigCmd = sprintf('%s\n\tfor v%i = 1:paramLengths(%i)', bigCmd, p, p);
    end
    bigCmd = sprintf('%s\n\t t=t+1;', bigCmd);
    for p=1:nParams
        bigCmd = sprintf('%s\n\t\t %s(t) = design.parameters.%s(v%i);', bigCmd, params{p}, params{p}, p);
    end
    for p=1:nParams
        bigCmd = sprintf('%s\nend', bigCmd);
    end
    
    try
        eval(bigCmd); 
    catch me
        keyboard
    end
    
    %then repmat nReps times 
    for p=1:nParams
        eval(sprintf('trials.%s = repmat(%s, numReps, 1);', params{p}, params{p}));
    end
    
    %shuffle order
    trials = trials(randperm(size(trials,1)),:);

    %cut off at the first numTrials
    trials=trials(1:numTrials,:);
    
end

%%%%%%%%%%%%%
% Make Random Variables
%%%%%%%%%%%%%

if isfield(design,'randVars')
    randParamNames=fieldnames(design.randVars);
    
    for p=1:length(randParamNames)
        pname=randParamNames{p};
        eval(sprintf('thesePs = design.randVars.%s;', pname));
        %make sure they're all column vectors
        if size(thesePs,1)==1
            thesePs = thesePs';
        end
        
        %how many possible values of this variable there are 
        nPs = length(thesePs); 
        
        %min. number of times each value will need to be repeated 
        nFullRep = floor(numTrials/nPs); 
        %beyond that, how many will need to be repeated one more time
        nExtra = numTrials - nFullRep*nPs; 
        
        %shuffle the order so that if there are some extra repetitions it's not always
        %the same ones
        shuffPs = thesePs(randperm(nPs));
        
        %create a vector with 1 row for each trial, using repmat and adding
        %on any extras
        theseVars = [repmat(thesePs, nFullRep, 1); shuffPs(1:nExtra)];
        
        %shuffle the order and add it to the trials structure 
        eval(sprintf('trials.%s = theseVars(randperm(numTrials));', pname));
        
  
    end
end

%% add trial number 
trials.trialNum = (1:numTrials)';





% thresh = estimateSIAM(ss)
%
% 2012 by Alex White
%
% estimate threshold from a SIAM staircase (Kaernbach, 1990)
% Input:
% - ss is the staircase structure
% - threshType is either 1 for average intensity levels on reversals
% (except for some specified reversals when step size was changing);
% or 2 for median of stimulus intensity levels on all trials except specified
% reversals.
%
% Output:
% - thresh is the threshold estimate
% - trialsCounted is the number of trials used in the estimate, after
% ignoring trials just after resets
%
% This takes the median of intensity levels used on trials in which the
% step size was stable.
%
% The staircase halves the step size after the first n reversals (ss.revsToHalveI), then
% keeps it stable until some number  reverals have happened
% (ss.resetRev), and then resets the counter and resets the step size, and
% then again halves the step size after each of the next n reversals.
% This function uses intensity values for all trials other than those for
% which the step size had just been reset or halved.

function [thresh, trialsCounted, threshSD] = estimateSIAM(ss, threshType)

%just take the median intensity value for all trials after we stopped
%halving the stimulus intensity level

%Figure out which trials to count

%Figure out which reversals to ignore
revToIgnore=[];
for rri = ss.resetRev
    revToIgnore = [revToIgnore rri rri+ss.revsToHalveI];
end
revToIgnore = revToIgnore(revToIgnore<=ss.nreversals);

if threshType == 1 %threshold is mean of intensities at reversals
    %this is unnecessary now that we keep track separately of intensities
    %at reversals with no change in step size
    %     revToCount = setdiff(1:ss.nreversals, revToIgnore);
    %     thresh = mean(ss.reversalInts(revToCount));
    %     threshSD = std(ss.reversalInts(revToCount));
    %     trialsCounted = length(revToCount);
    
    thresh = mean(ss.revStableStepsIs); 
    threshSD = std(ss.revStableStepsIs); 
    trialsCounted = length(ss.revStableStepsIs); 
    
elseif threshType == 2 %threshold is median of all intensities except after reversals when step size was changing
    
    %Ignore the trials between these ignored reversals
    tsToIgnore = [];
    for rti = revToIgnore
        if rti==0
            startt = 1;
        else
            startt = ss.reversalTs(rti);
        end
        if rti>=length(ss.reversalTs)
            endt = ss.tnum;
        else
            endt = ss.reversalTs(rti+1);
        end
        
        tsToIgnore = [tsToIgnore startt:endt];
    end
    
    tsToIgnore=unique(tsToIgnore);    
    tsToCount = setdiff(1:ss.tnum,tsToIgnore);
    
    thresh = median(ss.ints(tsToCount));
    threshSD = std(ss.ints(tsToCount));
    trialsCounted = length(tsToCount);
    
else
    fprintf(1,'\n\n(estimateSIAM) Error: threshType should be 1 (for reversal values) or 2 (for all intensity values)\n\n');
    return;
end


% ss = updateSIAM(ss,pres,resp)
% SIAM - single interval adjustment matrix procedure for threshold estimate
%
% 2012 by Alex White
%
% input:
% - ss, the SIAM staircase structure
% - pres, whether the signal was present (1) or absent (0) on this last trial
% - resp, whether the subjects said yes (1) or no (0)
%
% output:
% - ss, the updated SIAM staircase structure
%
%
% The staircase halves the step size after the first n reversals (ss.revsToHalveI), then
% keeps it stable until some number (ss.resetRev) reverals have happened, and then resets the
% counter and resets the step size, and then again halves the step size after each of the next 
% n reversals. 
% 
% It also resets the step size and counter if it seems to have gotten stuck - 
% that is, if the intensity values on the last 4 hits were the same. This
% can happen if the step size has become too small for the available
% intensity levels. 

function ss = updateSIAM(ss,pres,resp)

ss.pres = [ss.pres pres];
ss.resp = [ss.resp resp];
ss.ints = [ss.ints ss.intensity];
ss.tnum = ss.tnum+1;


if pres && resp         %hit
    step = ss.steps(1);
elseif pres && ~resp    %miss
    step = ss.steps(2);
elseif ~pres && resp    %false alarm
    step = ss.steps(3);
elseif ~pres && ~resp   %correct rejection
    step = ss.steps(4);
end

nextInt = ss.intensity+step;
% set intensity to one of the possible ones in the set
if ~isempty(ss.iSet)
    diffs = abs(nextInt-ss.iSet);
    nextInt = ss.iSet(diffs == min(diffs));
    if length(nextInt)>1
        %if there are two testimates because the desired level is between two
        %possible ones, then choose the one that's closer to the last value
        chgmag = abs(nextInt-ss.intensity);
        nextInt = nextInt(chgmag == min(chgmag));
    end
%clip intensity to set minimum and maximum
elseif ~isempty(ss.max)
    nextInt(nextInt<ss.min) = ss.min;
    nextInt(nextInt>ss.max) = ss.max;
end
    
ss.intensity = nextInt;

if step~=0
    dir = step/abs(step); %whether we're increasing or deacreasing intensity
else
    dir = ss.direction;
end

if ss.tnum>1
    reversal = dir ~= ss.direction;
    if reversal
        ss.reversalInts = [ss.reversalInts ss.ints(end)];
        ss.reversalTs = [ss.reversalTs ss.tnum];
        ss.nreversals = ss.nreversals+1;
        
        effectiveNRev = ss.nreversals - ss.resetRev(end);
        
        %Halve the step size?
        if any(effectiveNRev==ss.revsToHalveI)
            ss.steps=ss.steps/2;
        else %keep track of the number of "good" reversals on which step size was not changing 
            ss.nRevStableSteps = ss.nRevStableSteps+1;
            ss.revStableStepsTs = [ss.revStableStepsTs ss.tnum];
            ss.revStableStepsIs = [ss.revStableStepsIs ss.ints(end)];
        end
        
        %Reset the step size?
        if effectiveNRev==ss.revToReset
            ss.steps=ss.origSteps;
            
            %reset counter
            ss.resetT=[ss.resetT ss.tnum];
            ss.resetRev=[ss.resetRev ss.nreversals];
        end
        
    else %check if we've gotten stuck
        %if the last 4 hits have been at exactly the same level, reset
        %and call this a 'reversal,' because step is changing from 0 to
        %something else
        
        hitIs = find(ss.pres & ss.resp);
        if length(hitIs)>=ss.nStuckToReset && ss.resetT(end)<(ss.tnum-ss.nStuckToReset)
            if all(ss.ints(hitIs((end-ss.nStuckToReset+1):end)) == ss.intensity)
                
                ss.reversalInts = [ss.reversalInts ss.ints(end)];
                ss.reversalTs   = [ss.reversalTs ss.tnum];
                ss.nreversals   = ss.nreversals+1;
                
                ss.steps = ss.origSteps;
                
                %reset counter
                ss.resetT   = [ss.resetT ss.tnum];
                ss.resetRev = [ss.resetRev ss.nreversals];
            end
        end
    end
end
ss.direction=dir;



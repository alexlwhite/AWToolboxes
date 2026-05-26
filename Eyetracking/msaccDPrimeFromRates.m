%function [d] = msaccDPrimeFromRates(presRate, abstRate, dT, correctP)
%          by Alex White, April 2015
%
% Computes  microsaccade dprime given smoothed microsaccade rates on
% stimulus present and absetn trials, with the rates in Hz and computed
% every dT milliseconds. 
% 
% Inputs: 
% - presRate: microsaccade rate in Hz on target present trials
% - abstRate: microsaccade rate in Hz on target absent trials
% - dT: the separation between timepoints (in ms) at which the rate is computed. 
% - correctP: the minimum allowed difference a rate (HR or FAR) from 0 or 1. 
%   Explanation: If hit rate or false alarm rate is very near 1 (or 0), need to
%   subtract (or add) one 'hit' or 'false alarm' and assume some number of
%   trials, within some window of time within which the smoothing filter applies. 
%   Thid becomes an issue when msacc rate drops essentially to 0 for a time
%   period wider than the filter used to smooth the rate. Say W is that width 
%   (so far W has been ~200 ms). 
%   So if the msacc rate is near 0 for a while and the HR
%   goes to nearly 1, let's assume that out of N trials there was at
%   leats 1 microsaccade in that 200 ms window, and set the HR accordingly.
%   So, correctP should be set to 1/(W*N). Generally good to set N to twice
%   the number of target-present trials. 
% 
% Outputs: 
% - d: the smoothed dprime. 


% relation between hit rate (HR) and *target-present* msacc rate P 
% P is in Hz. So at time point t (1 ms), rate of P means that P/1000 msaccs
% happen every trial at time t. Or, if you ran 1000 trials, the total number
% of msaccs recorded at t would be P. In other words, the proprtion of trials with 
% msaccs at time t is P/1000 =   the miss rate. 
% Therefore, HR = 1-P/1000. 
% If P = 0, HR = 1. 

% Relation between false alarm rate (FAR) and *target-absent* msacc rate A:
% If you ran 1000 trials,the total number
% of msaccs recorded at t would be A. So, A/10 = correct reject rate. 
 %therefore, FAR = 1-A/1000. 
 
 %To avoid ceiling effects, lets say we add Q to both P and A. 
 %Then HR = 1-(P+Q)/1000; FAR = 1-(A+Q)/1000
 
 %Or if HR is too high, assum that 1 msacc happened in 1000 trials, which
 %would be P = 0.001 Hz, HR = 1-0.0001 = 0.9999


function [d, corrected] = msaccDPrimeFromRates(presRate, abstRate, dT, correctP)

%Assume rates are in Hz 
%Given rate in Hz, what's the proportion of trials with a msacc at each
%time t? 
pMsaccPres = presRate/(1000/dT); 
pMsaccAbst = abstRate/(1000/dT); 

%A "hit" is absence of microsaccade when target was present, so the hit rate 1
%minus the probability of a microsaccade on target present trials
hitR = 1-pMsaccPres;
%A "false alarm" is the absence of a microsaccade when the target was
%absent, for the false alarm rate is 1 minus the probability of msacc on
%target absent trials
FAR  = 1-pMsaccAbst;

%Avoid ceiling effects 
corrected = hitR>(1-correctP) | hitR<correctP | FAR>(1-correctP) | FAR<correctP;

hitR(hitR>(1-correctP)) = (1-correctP);
hitR(hitR<correctP) = correctP; 

FAR(FAR>(1-correctP)) = (1-correctP);
FAR(FAR<correctP) = correctP; 
    
d = norminv(hitR) - norminv(FAR);
            



    



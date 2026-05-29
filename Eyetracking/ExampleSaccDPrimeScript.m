
saccTimes = [ ]; %a vector of onset times of saccades, accumulated across trials of a given condition, relative to "time 0" when a stimulus or a blank was presented. 
numTrials = 100; %or whatever the number of trials over which these saccTimes were integrated 
rateTLims = [0 500]; % time window in which you want to compute the saccade rate, e.g., from 0 ms before to 700 ms after stimulus onset. 
kernelType = 3; %the type of smoothing kernel. I usually use type 3, the "causal" kernel 
kernelWidth = 20; %the width of the smoothing kernel in ms 
dT = 1; %number of ms between samples 

%It also helps to know how wide the filter is at the base, for telling the code how to correct hit or false alarm rates that hit 0 or 1
ts = 0:1:500;   t=250;
alpha = 1/kernelWidth;               
cts = t-ts+ip.kernelWidth; %adding KernelWidth shifts the filter back in time so its peak is at t
w = (alpha^2).*cts.*exp(-1*alpha*cts);
w(w<0) = 0;
maxH = max(w);
critH = maxH*0.1; %10% of total height cutoff 
filterBaseWidth = sum(w>=critH);

%  if "hit rate" gets to small, correct by assuming 1 microsaccade in a window filterBaseWidth ms wide, in twice the number of target present trials:                    
correctP =  1/(2*numTrials*filterBaseWidth);

%% Now compute the saccade rate for this condition 
[sR,rT] = saccRate(saccTimes,rateTLims,kernelType,kernelWidth,dT,numTrials);

%do that saccade rate computation  for both stimulus-present and stimulus absent (blank or baseline) conditions, to get "presRate" and "abstRate"

%% hen calculate oculomotor d' : 
saccDPrime = msaccDPrimeFromRates(presRate, abstRate, dT, correctP);

%% then take the maximum saccDPrime within a certain time window, and/or take the cumulative sum wtihin a certain time window. 


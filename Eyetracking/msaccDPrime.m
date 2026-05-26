% [d,dts,ratesCorrected] = msaccDPrime(onsetTimesPres,onsetTimesAbst,numTrialsPres,numTrialsAbst,tLims,kernelType,kernelWidth,dT)
%          by Alex White, March 2014
%
% Computes smoothed microsaccade dprime (or rate of any event), given a vector of onset
% times, with a particular filter kernel, within some specified time limits,
% and with a particular temporal resolution.
% 
% Inputs: 
% - onsetTimesPres: a vector of times at which event onsets occurred (in
% ms) when stimulus was present
% - onsetTimesAbst: a vector of times at which event onsets occurred (in
% ms) when stimulus was absent
% - numTrialsPres: the total number of trials included in onsetTimesPres.
% - numTrialsAbst: same, but for target absent trials 
% - tLims: the time limits [tMin tMax] over which to compute the smooted
%   rate. 
% - kernelType: 1=Boxcar of width kernelWidth; 2=Gaussian, with sigma=kernelWidth;
%               3=causal, with alpha = 1/kernelWidth
% - kernel width
% - dT: the separation between timepoints (in ms) at which the rate is computed. 
%   This is also width of the bins used for counting up the number of event onsets 
%   that occurred at each timepoint. 

% 
% Outputs: 
% - d: the smoothed dprime. The length of this vector is length(tMin:dT:tMax)
% - dts: a vector of times (in ms) for each element of r 
% - ratesCorrected: a vector with one element for each timepoint,
%   indicating whether or not either the hit rate or false alarm rate was
%   corrected to avoid undefined d'. 



function [d,dts,ratesCorrected] = msaccDPrime(onsetTimesPres,onsetTimesAbst,numTrialsPres,numTrialsAbst,tLims,kernelType,kernelWidth,dT)

%time limits for dprime computation. 
minT = tLims(1); % min([onsetTimes tLims(1)]); 
maxT = tLims(2); % max([onsetTimes tLims(2)]); 

%define specific times at which to compute dprime
ots = minT:dT:maxT; 
%bin edges: 
edges=(ots(1)-dT/2):dT:(ots(end)+dT/2); 
%bin centers: 
dts = mean([edges(1:(end-1)); edges(2:end)]); 

%count up the number of microsaccade onsets in each time bin with width dT, for
%target present trials 
countPres = histc(onsetTimesPres,edges); 
%exclude the last time bin which from histc is always 0, because it's the
%number of onsetTimes that match edges(end), which we set up to be beyond
%maxT. 
countPres = countPres(1:(end-1));
if isempty(countPres), countPres = zeros(size(dts)); end

%count up the humber of microsaccade onsets in each time bin, for target
%absent trials 
countAbst = histc(onsetTimesAbst,edges); 
%exclude the last time bin which from histc is always 0, because it's the
%number of onsetTimes that match edges(end), which we set up to be beyond
%maxT. 
countAbst = countAbst(1:(end-1));
if isempty(countAbst), countAbst = zeros(size(dts)); end


%for all time points t, compute hit and false alarm rates: 
%hit: no microsaccade when there was a stimulus
hitRs = 1-countPres/numTrialsPres;

%false alarm: lack of microsaccade when there was no stimulus
falRs = 1-countAbst/numTrialsAbst; 


smoothHitRs = zeros(size(dts));
smoothFalRs = zeros(size(dts));
ratesCorrected = false(size(dts));
ts = 1:length(dts);

if kernelType==3
    alpha = 1/kernelWidth;
end

%smooth (temporally average) hit and false alarm rates 
for t=ts
        
    %boxcar kernel 
    if kernelType==1
        w=zeros(1,length(ts)); 
        w(dts>=(dts(t)-kernelWidth/2) & dts<=(dts(t)+kernelWidth/2)) = 1; 
    %gaussian kernel     
    elseif kernelType==2
        w=normpdf(ts,t,kernelWidth);
    %'causal kernel'
    elseif kernelType==3
        cts = t-ts+kernelWidth; %adding KernelWidth shifts the filter back in time so its peak is at t
        w=(alpha^2).*cts.*exp(-1*alpha*cts);
        w(w<0)=0;
    end
    %normalize w to sum of 1
    w = w./sum(w); 
    
    smoothHitRs(t) = hitRs*w';
    smoothFalRs(t) = falRs*w';
    
end

%Avoid round-off errors 
smoothHitRs = round(smoothHitRs*10000000)/10000000;
smoothFalRs = round(smoothFalRs*10000000)/10000000;


%Deal with case when hit rate is 1 or 0

badIs = smoothHitRs==0; 
ratesCorrected(badIs) = true;
smoothHitRs(badIs) = 0.0001;

badIs = smoothHitRs==1; 
ratesCorrected(badIs) = true;
smoothHitRs(badIs) = 0.9999;

badIs = smoothFalRs==0; 
ratesCorrected(badIs) = true;
smoothFalRs(badIs) = 0.0001;

badIs = smoothFalRs==1; 
ratesCorrected(badIs) = true;
smoothFalRs(badIs) = 0.9999;

%don't allow dprime below zero
% badIs = smoothFalRs>smoothHitRs; 
% ratesCorrected(badIs) = true;
% smoothFalRs(badIs) = smoothHitRs(badIs); 
    
d = norminv(smoothHitRs) - norminv(smoothFalRs);
            


    



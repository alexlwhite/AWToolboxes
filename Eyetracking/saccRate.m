% [r,rT] = saccRate(onsetTimes,tLims,kernelType,kernelWidth,dT,numTrials)
%          by Alex White, March 2014
%
% Computes smoothed saccade rate (or rate of any event) as a function of time, 
% given a vector of onsettimes, with a particular filter kernel, within some specified time limits,
% and with a particular temporal resolution.
% 
% Inputs: 
% - onsetTimes: a vector of times at which event onsets occurred (in ms) 
% - tLims: the time limits [tMin tMax] over which to compute the smooted
%   rate. 
% - kernelType: 1=Boxcar of width kernelWidth; 2=Gaussian, with sigma=kernelWidth;
%               3=causal, with alpha = 1/kernelWidth
% - dT: the separation between timepoints (in ms) at which the rate is computed. 
%   This is also width of the bins used for counting up the number of event onsets 
%   that occurred at each timepoint. 
% - numTrials: the total number of trials included in onsetTimes 
% - doPlot (optional): whether or not to make a plot of rate computation
% with filter 
% 
% Outputs: 
% - r: the smoothed rate. The length of this vector is length(tMin:dT:tMax)
% - rT: a vector of times (in ms) for each element of r 



function [r,rT] = saccRate(onsetTimes,tLims,kernelType,kernelWidth,dT,numTrials,doPlot)

if ~exist('doPlot','var'), doPlot = false; end

%time limits for rate computation. 
minT = tLims(1);  
maxT = tLims(2); 

%count up the number of onsets in bins with width dT 
%first prepare bin edges
ots = minT:dT:maxT; 
%bin edges: 
edges=(ots(1)-dT/2):dT:(ots(end)+dT/2); 
%bin centers: 
rT = mean([edges(1:(end-1)); edges(2:end)]); 

count = histcounts(onsetTimes,edges); 
%because we now use histcounts and not histc, I don't think we need to do
%this: 
%exclude the last time bin which from histc is always 0, because it's the
%number of onsetTimes that match edges(end), which we set up to be beyond maxT. 
%count = count(1:(end-1));

if isempty(count), count = zeros(size(rT)); end

%smooth rate r
r = zeros(size(rT));
ts = 1:length(rT);

if kernelType==3
    alpha = 1/kernelWidth;
end

%for all time points t, compute dot-product between event count vector and kernel w 
for t=ts
    %boxcar kernel 
    if kernelType==1
        w = zeros(1,length(ts)); 
        w(rT>=(rT(t)-kernelWidth/2) & rT<=(rT(t)+kernelWidth/2)) = 1; 
    %gaussian kernel     
    elseif kernelType==2
        w = normpdf(ts,t,kernelWidth);
    %'causal kernel'
    elseif kernelType==3
        cts = t-ts+kernelWidth; %adding KernelWidth shifts the filter back in time so its peak is at t
        w = (alpha^2).*cts.*exp(-1*alpha*cts);
        w(w<0) = 0;
    end
    %normalize w to sum of 1
    w = w./sum(w); 
    
    %dot product 
    r(t) = count*w';

    %save w at timepoint closes to 0   
    if doPlot && t==ts(rT==min(abs(rT)))
        middleW = w; 
        middleTI = t;
    end
    
end

%divide by the number of trials
r = r./numTrials; 

%adjust for temporal resolution to get into units of Hz
r = (1000/dT)*r; 

%plot
if doPlot
    subplot(2,1,1); hold on;
    plot(rT,count,'g.','MarkerSize',8);
    plot(rT,middleW*max(count)/max(middleW),'b-'); %rescale w for visualization 
    plot(ones(1,2)*rT(middleTI),[0 max(count)],'k-');
    xlabel('Time rel to target onset');
    ylabel('Number of microsaccade onsets');
    ylim([0 max(count)+0.5]);
    
    subplot(2,1,2); hold on;
    plot(ots,r,'r-');
    xlabel('Time rel to target onset');
    ylabel('Microsaccade rate (Hz)'); 
end
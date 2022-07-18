%% function [medPos, meanPos, goodTimes, blinkCutTimes, noBlinkIntervals, pDataRemainAfterBlinkCut] = computeGazePosAndBlinks(time1, time2, edf)
% This function finds blinks in a sequence of gaze positions and computes
% the median gaze position (excluding periods with blinks). Blinks are
% detected because pupil size is 0.
%
% Inputs:
% - time1: starting time to analyze, in units used by edf.Samples.time
% - time2: ending time to analyze, also in units used by edf.Samples.time
% - edf: structure returned by Edf2Mat
%
%
% Outputs
% - medPos: 1x2 vector, median [horizontal, vertical] gaze position, in PIXELS (from edf.Samples.posX)
% - meanPos: 1x2 vecotr, [horizontal, vertical] gaze position, in PIXELS (from edf.Samples.posX)
% - goodTimes: a 1xT vector of Boolean values, which is true except during
% intervals with blinks. T is the number of samples in the interval to be
% analyzed.
% - blinkCutTimes: a Bx2 matrix that defines intervals to cut out because
%   of blnks. B is the number of blinks detected. The 1st column is start
%   times, and 2nd column is end times. These times are absolute values from
%   edf.Samples.time.
% - noBlinkIntervals: the converse of blinkCutTimes. A matrix of size B+1 x 2.
%   Each row is for a period with no blinks, [starttime endtime]. Again,
%   these times are absolute values from edf.Samples.time.
% - pDataRemainAfterBlinkCut: proportion of data remaining after cut out blinks.
%
% Note: the edf file aslo has blinksEvents.eBlink!
function [medPos, meanPos, goodTimes, blinkCutTimes, noBlinkIntervals, pDataRemainAfterBlinkCut] = computeGazePosAndBlinks(time1, time2, edf)

doPlot = false;

%some things are hard-coded: 
beforeBlinkBuffer = 100; %ms before blink starts to consider cutting (eventually depends on eye velocity)
afterBlinkBuffer = 120; %ms after blink starts to consider cutting (eventually depends on eye velocity)
maxYVelForCut = 0.5; % cutoff y-velocity in pix/ms to consider as when the eye starts or stops moving due to blink artifact 

intime = edf.Samples.time>=time1 & edf.Samples.time<time2;

times = edf.Samples.time(intime);
samrat = round(1000/mean(diff(times)));
msPerSample = 1000/samrat; %milliseconds per sample


eyeX = edf.Samples.posX(intime);
eyeY = edf.Samples.posY(intime);
pupSz = edf.Samples.pupilSize(intime);

%detect blinks as pupil size is 0 or NaN (different edf processing tools
%either insert a 0 or leave as NaN
isBlink = (pupSz <= 0) | isnan(pupSz);

%find onset and offset of blinks
blinkOn = false; blinkCount = 0;
blinkOnsetIs = []; blinkOffsetIs = [];
for tt = 1:length(times)
    if isBlink(tt)
        if ~blinkOn
            blinkOnsetIs = [blinkOnsetIs tt];
            blinkCount = blinkCount + 1;
        end
        blinkOn = true;
    else
        if blinkOn
            blinkOffsetIs = [blinkOffsetIs tt];
        end
        blinkOn = false;
    end
end
if blinkOn, blinkOffsetIs = [blinkOffsetIs tt]; end


goodTimes = true(size(times));

uniqueBlinkCount = 0;
blinkCutTimes = [];

minGapBtwnBlinks = 2; %gaps must be separated by this many ms, otherwise are merged
%determine which times are "good", cutting out the blink with some padding
%on either side:

if doPlot, newfig = figure; end

if blinkCount>0
    for bci=1:blinkCount
        
        minT = floor(blinkOnsetIs(bci)-beforeBlinkBuffer/msPerSample);
        minT(minT<1) = 1;
        
        maxT = ceil((blinkOffsetIs(bci)+afterBlinkBuffer/msPerSample));
        maxT(maxT>length(times)) = length(times);
        
        thisBlinkTimes = minT:maxT;
        
        yvels = [1; diff(eyeY(thisBlinkTimes))/msPerSample];
        zeroVelTimes = thisBlinkTimes(abs(yvels)<maxYVelForCut); %sometimes the velocity cross 0 but due to sampling doesnt quite reach 0, so we'll take 0.5 pix per ms as cutoff
        
        %time when the blink distortions started: latest time before signal
        %disappears when y velocity was 0
        minAroundShut = 10;
        zeroTimesBeforeShut = zeroVelTimes - blinkOnsetIs(bci);
        
        zeroTimesBeforeShut = zeroTimesBeforeShut(zeroTimesBeforeShut<(-minAroundShut));
        if isempty(zeroTimesBeforeShut)
            zeroTimesBeforeShut = min(zeroVelTimes - blinkOnsetIs(bci));
        end
        startCutBuffer = max(zeroTimesBeforeShut);
        if ~isempty(startCutBuffer)
            startCutT = blinkOnsetIs(bci)+startCutBuffer;
        else
            startCutT = blinkOnsetIs(bci);
        end
        
        %time when the blink distortions stopped: earliest time after
        %singal disappears when y velocity was 0
        zeroTimesAfterShut = zeroVelTimes - blinkOffsetIs(bci);
        zeroTimesAfterShut = zeroTimesAfterShut(zeroTimesAfterShut>(minAroundShut));
        if isempty(zeroTimesAfterShut)
            zeroTimesAfterShut = max(zeroVelTimes - blinkOffsetIs(bci));
        end
        
        endCutBuffer = min(zeroTimesAfterShut);
        if ~isempty(endCutBuffer)
            endCutT = blinkOffsetIs(bci)+endCutBuffer;
        else
            endCutT = blinkOffsetIs(bci);
        end
        
        if doPlot
            figure(newfig); clf;
            subplot(2,1,1); hold on;
            plot(thisBlinkTimes, eyeY(thisBlinkTimes), 'b-');
            plot([startCutT endCutT], eyeY([startCutT endCutT]), 'b.', 'MarkerSize', 8);
            xlim([minT maxT]);
            xlabel('time'); ylabel('y position');
            
            subplot(2,1,2); hold on;
            plot([minT maxT], [0 0], 'k-');
            plot(thisBlinkTimes, yvels, 'r-');
            plot([startCutT endCutT], yvels([find(thisBlinkTimes==startCutT) find(thisBlinkTimes==endCutT)]), 'r.',  'MarkerSize', 8);
            
            xlim([minT maxT]);
            xlabel('time'); ylabel('y velocity');
            pause(.25);
            
        end
        
        goodTimes(startCutT:endCutT) = false;
        
        %make sure to avoid having one blink cut time begin during a
        %previous one. Merge continguous blink intervals
        if uniqueBlinkCount>0
            %if the start of this blink interval is later than the ending
            %of the previous one, separated by a min gap, then add a new blnk interval
            if times(startCutT) > (blinkCutTimes(uniqueBlinkCount,2)+minGapBtwnBlinks)
                uniqueBlinkCount = uniqueBlinkCount+1;
                blinkCutTimes(uniqueBlinkCount, :) = [times(startCutT) times(endCutT)];
            else %otherwise, just set the end time of the previous blink interval to the end of this one, as they are continuous
                blinkCutTimes(uniqueBlinkCount, 2) = times(endCutT);
            end
        else
            uniqueBlinkCount = uniqueBlinkCount+1;
            try
                blinkCutTimes(uniqueBlinkCount, :) = [times(startCutT) times(endCutT)];
            catch
                keyboard
            end
        end
        
    end
end



%find sections of time between blinks to analyze
nBlinks = size(blinkCutTimes,1);
if nBlinks>0
    if blinkCutTimes(1,1)>time1
        noBlinkIntervals = [time1 blinkCutTimes(1,1)-1];
    else
        noBlinkIntervals = [];
    end
    
    for bi = 1:nBlinks
        if bi<nBlinks && blinkCutTimes(nBlinks,2)<time2
            noBlinkIntervals = [noBlinkIntervals; blinkCutTimes(bi,2)+1 blinkCutTimes(bi+1,1)-1];
        end
    end
    if blinkCutTimes(nBlinks,2)<(time2-1)
        noBlinkIntervals = [noBlinkIntervals; blinkCutTimes(nBlinks,2)+1 time2];
    end
else
    noBlinkIntervals = [time1 time2];
end

%exclude no blink intervals with 0 duration 
goodDurs = 1+diff(noBlinkIntervals,1,2);
if any(goodDurs<1)
    keyboard
end

pDataRemainAfterBlinkCut = mean(goodTimes);

%compute median gaze position, excluding blinks
medX = median(eyeX(goodTimes));
medY = median(eyeY(goodTimes));

medPos = [medX medY];

%and the mean
meanX = mean(eyeX(goodTimes));
meanY = mean(eyeY(goodTimes));

meanPos = [meanX meanY];



%% function [medPos, meanPos, goodTimes, blinkCutTimes, noBlinkIntervals, pDataRemainAfterBlinkCut] = computeBinocGazePosAndBlinks(time1, time2, edf)
% This function finds blinks in a sequence of gaze positions and computes
% the median gaze position (excluding periods with blinks). This version can take in either monocular or binocular gaze position data,
% so it estimates median gaze position for left and right eye.
% Blinks are detected because pupil size is 0. Or NaN if you want (set
% countNaNPupilAsBlink to true)
%
% Inputs:
% - time1: starting time to analyze, in units used by edf.Samples.time
% - time2: ending time to analyze, also in units used by edf.Samples.time
% - edf: structure returned by Edf2Mat
%
%
% Outputs
% - medPos: 2x2 vector, median gaze position, in PIXELS (from edf.Samples.posX).
%    Rows = [left eye; right eye]; Columns = [horizontal, vertical]
% - meanPos: 2x2 vector, mean gaze position, in PIXELS (from edf.Samples.posX)
%    Rows = [left eye; right eye]; Columns = [horizontal, vertical]
% - goodTimes: a 1xT vector of Boolean values, which is true except during
% intervals with blinks. T is the number of samples in the interval to be
% analyzed.
% - blinkCutTimes: a Bx2 matrix that defines intervals to cut out because
%   of blinks. B is the number of blinks detected. The 1st column is start
%   times, and 2nd column is end times. These times are absolute values from
%   edf.Samples.time.
% - noBlinkIntervals: the converse of blinkCutTimes. A matrix of size B+1 x 2.
%   Each row is for a period with no blinks, [starttime endtime]. Again,
%   these times are absolute values from edf.Samples.time.
% - pDataRemainAfterBlinkCut: proportion of data remaining after cut out blinks.
% - pupilMissingTimes: a Bx2 matrix that defines start and end times of each blink
% (when pupil size is missing.   These times are absolute values from edf.Samples.time.
% Note: the edf file aslo has blinksEvents.eBlink!
% - B: a table with 1 row per blink, with onset & offset times and cut buffers

 
function [medPos, meanPos, goodTimes, blinkCutTimes, noBlinkIntervals, pDataRemainAfterBlinkCut, pupilMissingTimes, B] = computeGazePosAndBlinks(time1, time2, edf)

doPlot = false;

%some things are hard-coded:
maxBeforeBlinkBuffer = 250; %ms before blink starts to consider cutting (eventually depends on eye velocity)
maxAfterBlinkBuffer = 250; %ms after blink starts to consider cutting (eventually depends on eye velocity)
minBeforeBuffer = 30; %min buffer time in ms to cut out before missing eye data 
minAfterBuffer = 70; %min buffer time in ms to cut out after missing eye data (longer than the before-blink buffer b/c it usually takes longer for the gaze to re-settle.)
minGapBtwnBlinks = 10; %gaps must be separated by this many ms, otherwise are merged
filtSD = 10; %SD of Gaussian smoothing kernel that is applied to the velocities, to detect when the gaze distortions start/stop around blinks 

%whether to count NaN pupil size as a blink, or just 0.
countNaNPupilAsBlink = false;

%maxYVelForCut = 0.25; % cutoff y-velocity in pix/ms to consider as when the eye starts or stops moving due to blink artifact

intime = edf.Samples.time>=time1 & edf.Samples.time<time2;

times = edf.Samples.time(intime);
samrat = round(1000/mean(diff(times)));
msPerSample = 1000/samrat; %milliseconds per sample

%convert max and min buffer times to samples
maxBeforeBlinkBuffer = round(maxBeforeBlinkBuffer/msPerSample);
maxAfterBlinkBuffer = round(maxAfterBlinkBuffer/msPerSample);
minBeforeBuffer = round(minBeforeBuffer/msPerSample);
minAfterBuffer = round(minAfterBuffer/msPerSample);
minGapBtwnBlinks = round(minGapBtwnBlinks/msPerSample);

eyeX = edf.Samples.posX(intime, :);
eyeY = edf.Samples.posY(intime, :);

isBinoc = size(eyeX, 2)>1;
%Compute gaze velocity with smoothing: 
if isBinoc
    eyes = 1:2;
    %here we are not using vecvel as it expects; it expects horiz and
    %vertical gaze position in 2 columns; here we're giving it left-eye and
    %right-eye vertcal gaze positions, but that's ok. 
    bothYVels = vecvel(eyeY, samrat, 2);

else
    eyes = 1;
    %cheat a bit with vecvel, give it the same y-position data twice 
    bothYVels = vecvel(repmat(eyeY, 1, 2), samrat, 2);
end

bothYVelsNoSmoth = bothYVels;

%Smooth the velocities: 
% MATLAB's smoothdata 'gaussian' window size is defined as roughly 6 * sigma 
% to capture 99.7% of the distribution.
filtSD_samp = round(filtSD/msPerSample);

window_size = 6 * filtSD_samp;

bothYVels = smoothdata(bothYVels, 'gaussian', window_size);
bothYVels(isnan(bothYVelsNoSmoth)) = NaN;   

%compute acceleration:  
bothYAcc = vecvel(bothYVels, samrat, 2);
yAccSD = std(bothYAcc,[],1,"omitmissing");
accCutuff = median(bothYAcc, 1, 'omitnan') + yAccSD;
accCutuff = median(abs(bothYAcc),'omitmissing') + yAccSD;
% plott = 1000;
% figure; 
% subplot(4,1,1); plot(times(1:plott), eyeY(1:plott,:)); title('Y Position');
% subplot(4,1,2); plot(times(1:plott), bothYVelsNoSmoth(1:plott,:));  title('Y Velocity');
% subplot(4,1,3); plot(times(1:plott), bothYVels(1:plott,:)); title('Smoothed Y Velocity');
% subplot(4,1,4), plot(times(1:plott), bothYAcc(1:plott,:)); title('Smoothed Y Acceleration');
% 
% keyboard

%pupil size 
pupSz = edf.Samples.pupilSize(intime, :);

%detect blinks as pupil size is 0 or NaN (different edf processing tools
%either insert a 0 or leave as NaN
isBlink = any(pupSz <= 0, 2);
if countNaNPupilAsBlink
    isBlink = isBlink | isnan(pupSz);
else
    if any(isnan(pupSz))
        keyboard
    end
end
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

pupilMissingTimes = [times(blinkOnsetIs) times(blinkOffsetIs)];


goodTimes = true(size(times));

uniqueBlinkCount = 0;
blinkCutTimes = [];



if doPlot & blinkCount>0, newfig = figure(3); clf; newfig2=figure(4); clf; end
%determine which times have "good" eye data, cutting out the blink with some padding
%on either side:

B = table;

if blinkCount>0
    for bci=1:blinkCount

        b = table;
        b.blinkNum = bci;
        
        %Define blink time with some buffer before and after:
        minT = floor(blinkOnsetIs(bci)-maxBeforeBlinkBuffer/msPerSample);
        minT(minT<1) = 1;

        maxT = ceil((blinkOffsetIs(bci)+maxAfterBlinkBuffer/msPerSample));
        maxT(maxT>length(times)) = length(times);

        thisBlinkTimes = minT:maxT;

        startCutTs = NaN(1,2);
        endCutTs = NaN(1,2);


        if isBinoc
            %check if both eyes blinked or if this was just missing data in
            %one eye
            thesePSz = pupSz(blinkOnsetIs(bci):blinkOffsetIs(bci),:);
            theseBlink = thesePSz==0;
            if countNaNPupilAsBlink
                theseBlink = theseBlink | isnan(thesePSz);
            end

            %realBlink: true if both eyes have pupil size missing at same
            %time. False if only one eye does.
            %If true, then we try to cut out a section around the blink,
            %based on distortions to vertical gaze position.
            %If false, we just cut out the section with missing data.
            realBlink = any(theseBlink(:,1) & theseBlink(:,2));

            if ~realBlink, fprintf(1,'\nBlink %i is not a "real" blink because not in both eyes\n', bci); end
        else
            realBlink = true;
        end


        b.realBlinkOrJustMissingData = realBlink;

        % if this seems to not be a real blink of both eyes, but rather a
        % period of missing data in 1 eye, then we just cutout the period
        % with that missing data.
        if ~realBlink
            startCutTs = blinkOnsetIs(bci);
            endCutTs = blinkOffsetIs(bci);
            %otherwise (including if we just have monocular data), we cut out
            %the missing data and a buffer period before and after to get rid
            %of artificats in the gaze positions:
        else
            for eye = eyes %loop thru eyes

              
                %% Set start time of blink-related distortions: 

                %Based on velocity: 
                %compute velocity of vertical gaze position
                %yvels = bothYVels(thisBlinkTimes, eye);

                %find when the y-velocity went to zero
                %zeroVelTimes = thisBlinkTimes(abs(yvels)<maxYVelForCut); %sometimes the velocity cross 0 but due to sampling doesnt quite reach 0, so we'll take 0.25 pix per ms as cutoff

                % latest time before pupil disappears when the y velocity was 0
                % zeroTimesBeforeShut = zeroVelTimes - blinkOnsetIs(bci);
                % 
                % zeroTimesBeforeShut = zeroTimesBeforeShut(zeroTimesBeforeShut<(-minBuffer));
                % if isempty(zeroTimesBeforeShut)
                %     zeroTimesBeforeShut = min(zeroVelTimes - blinkOnsetIs(bci));
                % end
                
                %based on acceleration
                 %find when the y-acceleration exceeded 0.5 SD
                yaccs = bothYAcc(thisBlinkTimes, eye);
                tooFastTimes = thisBlinkTimes(abs(yaccs)>accCutuff(eye));

                slowEnoughTimes = thisBlinkTimes(abs(yaccs)<accCutuff(eye));
                slowEnoughBeforeShut = slowEnoughTimes - blinkOnsetIs(bci);
                slowEnoughBeforeShut = slowEnoughBeforeShut(slowEnoughBeforeShut<(-minBeforeBuffer));
                if isempty(slowEnoughBeforeShut) %if none are less than the buffer, take the earlist one 
                    slowEnoughBeforeShut = min(slowEnoughBeforeShut);
                end
                %take the last time that the eye was accelerating slow enough 
                startCutBuffer = min([-minBeforeBuffer max(slowEnoughBeforeShut)]);

                %now set startCutT: the sample index number when we should
                %start cutting out data: 
                if ~isempty(startCutBuffer)
                    %but don't cut more than maxBeforeBlinkBuffer (take max of 2
                    %negative numbers)
                    startCutBuffer = max([startCutBuffer -maxBeforeBlinkBuffer]);
                    if startCutBuffer>0, keyboard; end

                    startCutT = blinkOnsetIs(bci)+startCutBuffer; %add a negative number
                    %can't cut out more data than we have:
                    if startCutT<1, startCutT =1; end
                else
                    startCutT = blinkOnsetIs(bci);
                end

                startCutTs(eye) = startCutT; %these are indices of samples 

                %% Set end time of blink-related distoartions: 
                %time when the blink distortions stopped: earliest time after
                %singal disappears when y velocity was 0
                % zeroTimesAfterShut = zeroVelTimes - blinkOffsetIs(bci);
                % zeroTimesAfterShut = zeroTimesAfterShut(zeroTimesAfterShut>(minBuffer));
                % if isempty(zeroTimesAfterShut)
                %     zeroTimesAfterShut = max(zeroVelTimes - blinkOffsetIs(bci));
                % end


                %OR based on acceleration
                % zeroTimesAfterShut = tooFastTimes - blinkOffsetIs(bci);
                % zeroTimesAfterShut = zeroTimesAfterShut(zeroTimesAfterShut>(minBuffer));
                % if isempty(zeroTimesAfterShut)
                %     zeroTimesAfterShut = max(tooFastTimes - blinkOffsetIs(bci));
                % end


                %based on first time eye slowed down
                zeroTimesAfterShut = slowEnoughTimes - blinkOffsetIs(bci);
                zeroTimesAfterShut = zeroTimesAfterShut(zeroTimesAfterShut>(minAfterBuffer));
                if isempty(zeroTimesAfterShut)
                    zeroTimesAfterShut = max(tooFastTimes - blinkOffsetIs(bci));
                end


                %take earliest time y-acceleration slowed won
                endCutBuffer = max([minAfterBuffer min(zeroTimesAfterShut)]); %must be positive
                if ~isempty(endCutBuffer)
                    %but don't cut more than maxAfterBlinkBuffer
                    endCutBuffer = min([endCutBuffer maxAfterBlinkBuffer]);

                    endCutT = blinkOffsetIs(bci)+endCutBuffer;
                    %can't cut out more data than we have:
                    if endCutT>size(eyeX,1), endCutT = length(times); end
                else
                    endCutT = blinkOffsetIs(bci);
                end
                endCutTs(eye) = endCutT;

            end
        end

        %cut out times that are bad for both eyes (earliest cut
        %start time and latest cut end time)
        startCutT = min(startCutTs);
        endCutT = max(endCutTs);
        

        if doPlot
            figure(newfig); clf;
            subplot(3,1,1); hold on;
            plot(thisBlinkTimes, eyeY(thisBlinkTimes, :), 'b-');
            plot([startCutT endCutT], eyeY([startCutT endCutT], :), 'b.', 'MarkerSize', 12);
            xlim([minT maxT]);
            xlabel('time'); ylabel('y position');

            subplot(3,1,2); hold on;
            plot([minT maxT], [0 0], 'k-');

            theseYVels = bothYVels(thisBlinkTimes,:);

            plot(thisBlinkTimes, theseYVels, 'r-');
            plot([startCutT endCutT], theseYVels([find(thisBlinkTimes==startCutT) find(thisBlinkTimes==endCutT)], :), 'r.',  'MarkerSize', 12);

            xlim([minT maxT]);
            xlabel('time'); ylabel('y velocity');
            

             subplot(3,1,3); hold on;
            plot([minT maxT], [0 0], 'k-');

             theseYAccs = bothYAcc(thisBlinkTimes,:);

            plot(thisBlinkTimes, theseYAccs, 'g-');
            plot([startCutT endCutT], theseYAccs([find(thisBlinkTimes==startCutT) find(thisBlinkTimes==endCutT)], :), 'g.',  'MarkerSize', 12);
            for eyei=eyes
                plot([minT maxT], accCutuff([eyei eyei]), 'g:');
            end
            xlim([minT maxT]);
            xlabel('time'); ylabel('y acceleration');

            pause; %(0.5);

        end

        goodTimes(startCutT:endCutT) = false;

        %set blinkCutTimes in units of ms, using the edf file's timestamps.
        %blinkCutTimes is a B x 2 matrix for B blinks, start times and end
        %times
       
        %make sure to avoid having one blink cut time begin during a
        %previous one. Merge continguous blink intervals
        if uniqueBlinkCount>0
            %if the start of this blink interval is later than the ending
            %of the previous one, separated by a min gap, then add a new blnk interval
            if times(startCutT) > (blinkCutTimes(uniqueBlinkCount,2)+minGapBtwnBlinks)
                uniqueBlinkCount = uniqueBlinkCount+1;
                blinkCutTimes(uniqueBlinkCount, :) = [times(startCutT) times(endCutT)];
                addToTable = true;
            else %otherwise, just set the end time of the previous blink interval to the end of this one, as they are continuous
                blinkCutTimes(uniqueBlinkCount, 2) = times(endCutT);
                addToTable = false;
            end
        else %if this is the first one
            uniqueBlinkCount = uniqueBlinkCount+1;
            blinkCutTimes(uniqueBlinkCount, :) = [times(startCutT) times(endCutT)];
            addToTable = true;
        end

        if addToTable
            b.onsetTime = times(blinkOnsetIs(bci));
            b.offsetTime = times(blinkOffsetIs(bci));
            b.preBlinkBufferCut = blinkOnsetIs(bci) - startCutT;
            b.postBlinkBufferCut = endCutT - blinkOffsetIs(bci);

            B = [B; b];
        end
      
    end
end



%% find sections of time between blinks to that are safe to analyze 
nBlinks = size(blinkCutTimes,1);
if nBlinks>0
    %start onset of 1st good interval: 
    %if the first blink ocurred later than onset of the whole interval,
    %then the 1st safe interval is from time1 to start of 1st blink
    if blinkCutTimes(1,1)>time1
        noBlinkIntervals = [time1 blinkCutTimes(1,1)-1];
    else %otherwise we have to wait until after the 1st blink 
        noBlinkIntervals = [];
    end

    %now loop through blinks and add the interveing times to
    %noBlinkIntervals
    for bi = 1:(nBlinks-1) %up until the last blink 
        noBlinkIntervals = [noBlinkIntervals; blinkCutTimes(bi,2)+1 blinkCutTimes(bi+1,1)-1];
    end
    %now the time after the last blink: 
    if blinkCutTimes(nBlinks,2)<(time2-1)
        noBlinkIntervals = [noBlinkIntervals; blinkCutTimes(nBlinks,2)+1 time2];
    end
else %if there are no blinks, then the whole interval given to this function is safe to analzye 
    noBlinkIntervals = [time1 time2];
end

%exclude no blink intervals with 0 duration
goodDurs = 1+diff(noBlinkIntervals,1,2);
if any(goodDurs<1)
    keyboard
end

pDataRemainAfterBlinkCut = mean(goodTimes);

%compute median gaze position, excluding blinks
medX = median(eyeX(goodTimes, :), 1);
medY = median(eyeY(goodTimes, :), 1);

medPos = [medX' medY']; %%    Rows = [left eye; right eye]; Columns = [horizontal, vertical]

%and the mean gaze position
meanX = mean(eyeX(goodTimes, :),1);
meanY = mean(eyeY(goodTimes, :), 1);

meanPos = [meanX' meanY']; %Rows = [left eye; right eye]; Columns = [horizontal, vertical]

if any(isnan(medPos)) | any(isnan(meanPos))
    keyboard
end
%% plot
if doPlot & nBlinks>0
    figure(newfig2); hold on;
    %plot(times, eyeX, 'b-');
    plot(times, eyeY, 'b-');
    plot(times, pupSz, 'g-');

    ylims = [0 max([max(pupSz) max(eyeY(:))])];

    for bi=1:size(pupilMissingTimes,1)
        plot(pupilMissingTimes([bi bi],1), ylims, 'k-');
        plot(pupilMissingTimes([bi bi],2), ylims, 'k-');
        plot(pupilMissingTimes(bi,:), ylims([2 2]), 'k-');
    end
    for bi=1:nBlinks
        plot(blinkCutTimes([bi bi],1), ylims, 'r-');
        plot(blinkCutTimes([bi bi],2), ylims, 'r-');
        plot(blinkCutTimes(bi,:), ylims([2 2]), 'r-');
    end

    ylims = ylims*1.09;
    if diff(ylims)>0
        ylim(ylims);
    end

    keyboard
end
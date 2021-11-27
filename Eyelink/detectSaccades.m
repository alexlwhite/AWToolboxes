function saccTable = detectSaccades(postns,vel,velThresh,MINDUR,mergeInterval)
%-------------------------------------------------------------------
%
%  Adapted by Alex White from the function "microsacc.m", (Version 2.1, 03 OCT 05)
%  Detection of monocular candidates for microsaccades;
%  Please cite: Engbert, R., & Mergenthaler, K. (2006) Microsaccades
%  are triggered by low retinal image slip. Proceedings of the National
%  Academy of Sciences of the United States of America, 103: 7192-7197.
%
%
%-------------------------------------------------------------------
%
%  INPUTS:
%
%  postns:        a Nx2 vector of gaze positions, 1st column for horiz and  second for vert. Units could be deg or pix
%  vel:           a Nx2 vector of velocities, 1st column for horiz and second for vert
%  velThresh      a 1x2 vector of horizontal and vertical velocity thresholds
%  MINDUR           minimal saccade duration, in *samples*
%  mergeInterval: for subsequent saccade candidates, in *samples*
%
%  OUTPUTs:
%  saccTable, a table with 1 row per saccade. Colums:
%    - onsetSample, time of saccade onset in samples;
%    - offsetSample, time of saccade offset in samples;
%    - peakVelocity, in deg/s (or pix/s, if input postns is in pixels)
%    - startX, starting horizontal position (in same coordinate frame as input postns)
%    - endX, ending horizontal position
%    - startY, starting vertical position
%    - endY, ending vertical position
%    - dx, horizontal difference in start and end points
%    - dy, vertical difference in start and end points
%    - totalAmpX, maximum horizontal difference in gaze positions during the whole high-velocity event
%    - totalAmpY, maximum vertical difference during the whole event
%    - amp: amplitude of saccade, distance from starting to ending point 
%    - curveRatio: max deviation from straight line, devided by amp

%---------------------------------------------------------------------

doPlot = false;

%now, define an threshold ellipse in 2D velocity space. That is set by
%input parameter velThresh
radiusX = velThresh(1);
radiusY = velThresh(2);

% compute test criterion: ellipse equation:
%test>1 if the velocity exceeds the criterion
test = (vel(:,1)/radiusX).^2 + (vel(:,2)/radiusY).^2;
%indx: times when the velocity exceeded the criterion
indx = find(test>1);

% determine saccades
N = length(indx); %number of samples that had velocity over threshold
sac = [];
nsac = 0;
dur = 1; %counter of saccade duration (samples)
a = 1; %index of saccade start time
k = 1; %index of which sample we're at
while k<N
    if indx(k+1)-indx(k)==1 %if there is no gap in time
        dur = dur + 1; %then the current "event's" duration increases
    else %we've jumped time, so the previous period of high velocity ended
        if dur>=MINDUR %if this high-velocity event was long enough
            nsac = nsac + 1;
            b = k; %saccade offset time
            sac(nsac,:) = [indx(a) indx(b)]; %[starttime, endtime] (samples)
        end
        a = k+1; %move on
        dur = 1; %restart counter
    end
    k = k + 1;
end

% check for one last saccade at the very end, for minimum duration
if dur>=MINDUR
    nsac = nsac + 1;
    b = k;
    sac(nsac,:) = [indx(a) indx(b)];
end

% merge saccades that are separated by not enough time
if ~isempty(sac)
    sacc = sac(1,:);    % merged saccade matrix;; starting with 1st saccade
    s    = 1;           % index of saccades in sac
    sss  = 1;           % boolean for still same saccade
    nsac = 1;           % number of saccades after merge
    while s<size(sac,1)
        if ~sss
            nsac = nsac + 1;
            sacc(nsac,:) = sac(s,:);
        end
        if sac(s+1,1)-sac(s,2) <= mergeInterval
            sacc(nsac,2) = sac(s+1,2);
            sss = 1;
        else
            sss = 0;
        end
        s = s+1;
    end
    if ~sss
        nsac = nsac + 1;
        sacc(nsac,:) = sac(s,:);
    end
else
    sacc = [];
    nsac = 0;
end

%% convert all this into a table  (added by ALW, Feb 2020)
saccTable = table;

if nsac>0
    saccTable.onsetSample = sacc(:,1);
    saccTable.offsetSample = sacc(:,2);
    
    emptyVec = NaN(size(saccTable.onsetSample));
    saccTable.peakVelocity  = emptyVec;
    saccTable.startX        = emptyVec;
    saccTable.endX          = emptyVec;
    saccTable.startY        = emptyVec;
    saccTable.endY          = emptyVec;
    saccTable.dx            = emptyVec;
    saccTable.dy            = emptyVec;
    saccTable.totalAmpX     = emptyVec;
    saccTable.totalAmpY     = emptyVec;
    
    % compute peak velocity, horizonal and vertical components
    for s=1:nsac
        % onset and offset
        a = sacc(s,1);
        b = sacc(s,2);
        % saccade peak velocity (vpeak)
        vpeak = max( sqrt( vel(a:b,1).^2 + vel(a:b,2).^2 ) );
        
        saccTable.peakVelocity(s) = vpeak;
        
        %start and end points
        saccTable.startX(s) = postns(a,1);
        saccTable.startY(s) = postns(a,2);
        saccTable.endX(s)   = postns(b,1);
        saccTable.endY(s)   = postns(b,2);
        
        % saccade vector (dx,dy). These are the differences in x and y positions
        % from saccade start and end times (defined by velocity)
        saccTable.dx(s) = postns(b,1)-postns(a,1);
        saccTable.dy(s) = postns(b,2)-postns(a,2);
        
        % saccade amplitude (totalAmpX,totalAmpY): These are the
        % differences in minimum and maximum x and y positions across the
        % entire duration of the saccade
        is = sacc(s,1):sacc(s,2);
        [minx, ix1] = min(postns(is,1));
        [maxx, ix2] = max(postns(is,1));
        [miny, iy1] = min(postns(is,2));
        [maxy, iy2] = max(postns(is,2));
        dX = sign(ix2-ix1)*(maxx-minx);
        dY = sign(iy2-iy1)*(maxy-miny);
        
        saccTable.totalAmpX(s) = dX;
        saccTable.totalAmpY(s) = dY;
        
        % saccade curvature: at each point, compute distance of gaze
        % position from the straight line connecting start and end
        % positions. Then blur those a bit, and take the max distance 
        slope = saccTable.dy(s)/saccTable.dx(s);
        if isinf(slope)
            slope = saccTable.dy(s)/0.1;
        end
        intercept = saccTable.startY(s)-saccTable.startX(s)*slope;
        
        curveDists = distanceFromPointToLine(postns(is,1),postns(is,2),slope,intercept);
        
        %smooth to reduce noise 
        filtWid=4;
        filt = ones(1,filtWid)/filtWid;
        smoothDists = conv(curveDists, filt, 'same');
        
        saccTable.maxCurveDeviation(s) = max(smoothDists);
        
        amp = sqrt(saccTable.dx(s).^2 + saccTable.dy(s).^2);
        curveRatio = saccTable.maxCurveDeviation(s)/amp; 
        
        %thresh = 0.15;
        %doPlot = curveRatio>thresh;
        
        
        if doPlot
            figure; subplot(2,1,1); hold on;
            xs=[saccTable.startX(s) saccTable.endX(s)];
            ys=[saccTable.startY(s) saccTable.endY(s)];
            predys = xs*slope+intercept;
            plot(xs, ys, 'bo');
            plot(xs,predys,'k-');
            plot(postns(is,1), postns(is,2), 'r.');
            
            subplot(2,1,2); hold on;
            plot(is, curveDists, 'g.-');
            plot(is, smoothDists, 'r.-');
            keyboard
        end
       
        
    end
    
    %amp: amplitude of saccade, distance from starting to ending point 
    saccTable.amp = sqrt(saccTable.dx.^2 + saccTable.dy.^2);
    
    
    %curve ratio: max deviation from straight line, devided by amp
    saccTable.curveRatio = saccTable.maxCurveDeviation./saccTable.amp;

    
end

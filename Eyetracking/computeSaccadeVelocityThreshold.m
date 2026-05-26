%% function velThresh = computeSaccadeVelocityThreshold(vel, velThreshSDs)
% This function computes the velocity threshold for saccade detection,
% given all the velocities observed ("vel") and a single number,
% velThreshSDs, which is how many median-based velocity SDs the threshold
% should be. 
% 
% Inputs: 
% - vel: a Tx2 matrix of gaze position velocities, where the 1st column is horizontal and the 2nd column is vertical 
% - velThreshSDs: a single number. The threshold velocity will be the
%    median-based standard deviation times velThreshSDs. 
% 
% Output: 
% - velThresh: a 1x2 vector of [horizontal vertical] velocity thresholds. 
%   This effectively defines an ellipse in 2D velocity space. A saccade can
%   then defined as a deviation outside that ellipse, for a given duration.
%   
% 
function velThresh = computeSaccadeVelocityThreshold(vel, velThreshSDs)

%first, the median-based standard deviations of velocities in x- and y-directions:
%msdx = sqrt(   median(velX^2) - median(velX)^2  )

msdx = sqrt( median(vel(:,1).^2) - (median(vel(:,1)))^2 );
msdy = sqrt( median(vel(:,2).^2) - (median(vel(:,2)))^2 );

%if it's too small, use the mean rather than median 
if msdx<realmin
    msdx = sqrt( mean(vel(:,1).^2) - (mean(vel(:,1)))^2 );
    %if msdx<realmin  %% CAUTION: I (Martin?) took this out because the only critical
    %                 %% time when this can occur is for full-time blinks. I take these out
    %                 %% later in the analysis.
    %    error('msdx<realmin in microsacc.m');
    %end
end
if msdy<realmin
    msdy = sqrt( mean(vel(:,2).^2) - (mean(vel(:,2)))^2 );
    %if msdy<realmin  %% CAUTION: I took this out because the only critical
    %                 %% time when this can occur is for full-time blinks. I take these out
    %                 %% later in the analysis.
    %    error('msdy<realmin in microsacc.m');
    %end
end

%now, define an threshold ellipse in 2D velocity space: median SD times
%some factor, velThreshSDs
threshX = velThreshSDs*msdx;
threshY = velThreshSDs*msdy;

velThresh = [threshX threshY];
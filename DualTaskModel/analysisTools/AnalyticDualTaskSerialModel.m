%% function [pProcessBoth, pTask1First, slope, intercept, distFromAllOrNone, spokeDistFromAllOrNone] = AnalyticDualTaskSerialModel(singleAccs, dualAccs, doPlot)
% by Alex White, 2017
%
% Given accuracy (p(correct) or area under ROC curve) in dual and single-task 
% conditions for two tasks (e.g., left and right) this model finds the best-fitting
% generalized serial model. This is an analytic solution, by which I mean
% that it does not do any optimization search or simulate individual
% trials. It basically just re-parameterizes the dual-task accuracy levels
% into two other parameters, given single-task accuracy levlels. 
% 
% The model assumes that on each trial one task (or side) is attended 'first', and the 
% other second. On some fraction of trials, both get processed, and on the 
% remaining trials, only the first one is processed and the subject must guess about the other. 
% The proportion of trials when both are processed is called pProcessBoth.
% The proportion of trials when task number 1 is attended first is called
% pTask1First.
% 
%
% This function returns those two free parameters as well as a slope and
% intercept for the best-fitting line in the AOC plot. That line is
% determined completely by the data and the processBoth parameter. 
% 
% Inputs: 
% - singleAccs: a 1x2 vector of accuracy in single-task conditions for the
%   two tasks. Units: p(correct) or area under ROC curve. 
%   The 1st element goes on the y axis, and the second goes on the X axis.  
% - dualAccs: a 1x2 vector of accuracy in dual-task conditions for the same
%   two tasks. The 1st element goes on the y axis, and the second goes on the X axis.  
% - doPlot: whether to make an ROC plot showing data and model predictions.
% 
% Outputs: 
% - pProcessBoth: the proportion of trials when the subject successfully processed both 
%   stimuli [each with   p(correct) = single-task p(correct)]. On the remaining 
%   trials, the subject is forced completely guess about the task not attended 'first'. 
% - pTask1First: best-fitting parameter for proportion of dual-task trials when task number 1
%   is attended 'first'. 
% - slope: slope of the best-fitting line in the AOC plot
% - intercept: y-intercept of the best-fitting line in the AOC plot 
% - distFromAllOrNone: distance between dual-task data point and the
%   nearest point on the all-or-none switching model. Negative if data point
%   is below the line. 
% - spokeDistFromAllOrNone: a normalized measure of distance from serial
%   model, along the line that connects the dual-task accuracy point to the
%   unlimited capacity parallel model. This distance is normalized by (i.e., divided by) 
%   by the total length of the segement of that line betwen the unlimited
%   capacity point and the serial model prediction's line. This measure is
%   negative if accuracy is worse than the serial model predicts. This
%   measure can only be computed if both singleAccs<dualAccs. If not, it is
%   NaN. 
%
function [pProcessBoth, pTask1First, slope, intercept, distFromAllOrNone, spokeDistFromAllOrNone] = AnalyticDualTaskSerialModel(singleAccs, dualAccs, doPlot, chanceLevel)

if nargin<4
    chanceLevel = 0.5; 
end

%To calculate things like slopes and intercepts and
%distances in the AOC space, we need to subtract chance level (usually 0.5) from the accuracy levels and
%pretend like the axis limits really are [0 0.5]. 
singleAccs = singleAccs - chanceLevel; 
dualAccs = dualAccs - chanceLevel; 


%pull out data:
dualY  = dualAccs(1); %y-value
dualX = dualAccs(2); %x-value
singleY  = singleAccs(1);%y-value
singleX = singleAccs(2);%x-value


%Given single-task accuracy on the two tasks, singleY and singleX,
%and dual-task accuries, dualY and dualX,
%we want to solve for the dual-task model parameters pProcessBoth and
%pTask1First.

%First we can solve for the slope of all the dual-task lines
%slope = (.5-singleY)/(singleX-0.5); %this was true before we subtracted 0.5 from all data points 
slope = -singleY/singleX;

%Second, we can solve for the y-intercept of the best-fitting line, using our measured
%left and right dual-task performance
intercept = dualY - slope*dualX;

%Then we can find the predicted value of left-dual task performance when attention is
%devoted entirely to task 2, and task 1 is supposed to be ignored on 100%
%of trials. Call this dual1Ignored (dual-task one when ignored). It is the y-value of the line when the x-value is
%singleX.
dual1Ignored = slope*singleX+intercept;

% from this we can solve for the proportion of trials when only 1 side is
% processed and the other must be guessed. 
% We assume that accuracy when guessing is 0.5.
guessingPCorr = 0; %0.5; %no longer 0.5 because we subtracted that out 
pProcessOnlyOne = (dual1Ignored - singleY)/(guessingPCorr - singleY);
%clip this propability to be between 0 and 1
pProcessOnlyOne(pProcessOnlyOne<0)=0; 
pProcessOnlyOne(pProcessOnlyOne>1)=1; 

pProcessBoth  = 1-pProcessOnlyOne;

%[We could also do this using the predicted right-side accuracy when right side is
%supposed to be ignored, dual2Ignored
%dual2Ignored = (singleY-intercept)/slope;
%pProcessOnlyOne_B = (dual2Ignored - singleX)/(guessingPCorr - singleX); %]

% Then solve for pAttendTask1, or pAL (proportion of dual-task trials
% attend to left first)
pTask1First = (dualX - singleX)/(pProcessOnlyOne*guessingPCorr - singleX*pProcessOnlyOne);

%clip in case dual-task performance exceeds any serial model 
pTask1First(pTask1First>1) = 1;
pTask1First(pTask1First<0) = 0; 


%% find the minimum distance between the dual-task data point and the all-or-none serial prediction line 
%difference in y-intercepts of best-fitting line and all-or-none line
% dY = intercept - singleY;
% 
% %compute x-intercept of best-fitting line
% xIntercept = -1*intercept/slope;
% %difference between this and the all-or-none's x-intercept (which is simple
% %single-task accuracy level on x-axis)
% dX = xIntercept - singleX;
% 
% %hypotenuese
% h = sqrt(dX^2 + dY^2); 
% %some trig:
% distFromAllOrNone1 = sqrt(dX^2 - (h/2)^2);

distFromAllOrNone = distanceFromPointToLine(dualX,dualY,slope,singleY);

%if performance is WORSE than all-or-none model, make this negative
%difference in y-intercepts of best-fitting line and all-or-none line
dY = intercept - singleY;
if dY<0
    distFromAllOrNone = distFromAllOrNone*-1;
end

%% compute the "spoke distance": another measure of how serial vs parallel this performance level is: 
% how far the dual-task point is from the serial model prediction, along a
% line that connects the dual-task point to the unlimited capacity point 
% This only works if dual-task accuracy is worse on both sides than
% single-task accuracy. 

canDoNormDist = dualY<singleY && dualX<singleX;

if canDoNormDist

    %define that line
    slope2 = (singleY-dualY)/(singleX-dualX);
    intercept2 = dualY - slope2*dualX; %I = y-slope*x;

    intercept1 = singleY; %y-intersect of the all-or-none serial model line

    %find where that line intersects the serial model
    intersectX = (intercept2 - intercept1)/(slope-slope2);

    intersectY = slope2*intersectX + intercept2;
    intersectY2 = slope*intersectX + intercept1;
    if abs(intersectY-intersectY2)>(10^-13)
        keyboard
    end

    %now find the distance of dual-task data point from that point on the
    %serial model line (the one that is on the same line that ocnnects the data
    %point to the unlimited capacity parallel point)
    dist2 = sqrt((dualY - intersectY)^2 + (dualX - intersectX)^2);
    %if performance is WORSE than all-or-none model, make this negative
    if dY<0
        dist2 = -1*dist2;
    end
    %and distance from that point to unlimited capacity point
    wholeDist =  sqrt((singleY - intersectY)^2 + (singleX - intersectX)^2);

    %normalize distance from serial model by that total distance to the
    %unlimited capacity model
    spokeDistFromAllOrNone = dist2/wholeDist;
elseif all(dualAccs==singleAccs)
    spokeDistFromAllOrNone = 1;
else
    spokeDistFromAllOrNone = NaN;
end

%% plot
if doPlot
    dataMarkSz = 30;
    modelMarkSz = 10;
    
    figure; hold on;
    
    %add 0.5 back to everything 
    singleY = singleY + chanceLevel;
    singleX = singleX + chanceLevel;
    dualY = dualY + chanceLevel;
    dualX = dualX + chanceLevel;
    
    %plot box constrained by indepenent processing:
    plot([singleX singleX],[chanceLevel singleY],'k-');
    plot([chanceLevel singleX],[singleY singleY],'k-');
    
    %plot all-or-none serial model. Solve for that equation using x,y
    %coordinates of left single task accuracy:
    x0 = chanceLevel; y0=singleY;
    betaAllOrNone = y0 - slope*chanceLevel;
    xs = [.5 singleX];
    ys = slope*xs+betaAllOrNone;
    plot(xs,ys,'k-');
    
    %plot the best-fitting serial model line (fit with pProcessOnlyOne, still free
    %parameter of pAttendTask1)
    dual2Ignored = (singleY-intercept)/slope;

    xs = [dual2Ignored singleX];
    ys = slope*xs+intercept;
    plot(xs+chanceLevel,ys+chanceLevel,'r-');
    
    
    %Plot prediction of fixed-capacity parallel processing
    singleAcc = [singleX singleY];
    singleDs = sqrt(2)*norminv(singleAcc); %convert to d' (assuming A' is like PC?')
    pSamplesOnTask1 = 0:0.01:1;  %vector of proportion of fixed number of sensory "samples" devited to task 1
    
    % How to calculate d' for diffent numbers of "samples"
    % given sample sizes ss1 and ss2: ss1 = ss2/2
    % std1 = sqrt(2)*std2
    % generalized:
    %if ss1 = ss2*q, then
    %std1 = sqrt(1/q)*std2
    
    %therefore, d' relation is:
    % d1 = d2/sqrt(1/q)
    
    %dual-task dprimes
    fixedD1s = singleDs(1)./sqrt(1./pSamplesOnTask1);
    fixedD2s = singleDs(2)./sqrt(1./(1-pSamplesOnTask1));
    
    %dual-task AROCs,
    fixedA1s = normcdf(fixedD1s/sqrt(2));
    fixedA2s = normcdf(fixedD2s/sqrt(2));
    
    plot(fixedA1s,fixedA2s,'y-','LineWidth',2)
    
    
    %plot data:
    plot(chanceLevel,singleY,'b.-','MarkerSize',dataMarkSz);
    plot(singleX,chanceLevel,'b.-','MarkerSize',dataMarkSz);
    plot(dualX,dualY,'k.','MarkerSize',dataMarkSz)
    
    
    %plot the predicted dual-task performance given our 2 parameters:
    predDual = SerialDualTaskAccGivenSingleTask([singleY singleX],pTask1First,pProcessBoth);
    plot(predDual(2), predDual(1),'r.','MarkerSize',modelMarkSz);

    %plot the line that connects dual-task point to unlimited capacity 
    if canDoNormDist
        someXs = [intersectX singleX];
        someYs = slope2*someXs + intercept2;
        plot(someXs+chanceLevel, someYs+chanceLevel, 'g--');
    end
    axis square;
    xlabel('Right side p(c)');
    ylabel('Left side p(c)');
    xlim([chanceLevel 1]); ylim([chanceLevel 1]);
    ticks = chanceLevel:0.1:1;
    set(gca,'XTick',ticks,'YTick',ticks);
    
end
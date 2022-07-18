%function [pProcessBoth, pTask1First, hSing] =  plotAOCWithPredictions(as,es,plotOpt)
%
% This function plots an AOC based on single- and dual-task accuracy. It
% also fits a "generalized" serial model to the data. The diagonal line in
% the AOC is the "all-or-none" serial model which assumes that only 1
% stimulus is processed per trial and the observer gets no information at
% all about the other. The "generalized" serial model assumes that on some
% variable propotion of dual-task trials, both sides are processed (to the same degree
% as in the single-task conditions). This parameter (pProcessBoth) can be
% estimated from the input data. The serial model can also estimate the
% proportion of dual-task trials in which stimulus 1 is processed "first."
% The serial model assumes that one stimulus is processed first, and then
% on pProcessBoth proportion of trials, the other stimulus is processed as well. 
% This parameter, pTask1First, can account for the dual-task point not
% falling directly in the middle of the diagonal line, but shifted to one
% side or the other. 
% 
%Inputs:
%- as: 2x2 matrix of area-under-ROC-curve measures.
%    rows = single-task, dual-task
%    columns = left side, right side
%- es: a 2x2 matrix of error bars for each of the conditions in as
%- plotOpt: structure with fields including::
%      - edgeColors: a 2x3 matrix with each row a color for marker edge.
%                    rows = single-task, dual-task
%      - fillColors: like edgeColors, but for fill of markers
%      - markSz: size of markers
%      - axlims: axis limits (for both x and y). Default: [0.5 1]
%      - doLegend: whether to make legend
%      - doAxLabels: whether to add x and y-axis labels
%      - units: what the axis units are (AUC, p(correct), etc))
%      - plotSerialPrediction: boolean, whether to plot the best-fitting serial-model's
%        prediction for dual-task accurayc as a red point 
%      - sideLabels: a 1x2 cell array for the labels of the two axes
%      - chanceLevel: what chance level is
%
%
% Outputs: 
% - pProcessBoth: parameter of generalized serial model: the proportion of
%   trials in which both stimuli are processed 
% - pTask1First: parameter of generalized serial model: proportion of
%   trials when task 1 is processed 'first', meaning that only task 1 is
%   procesed on (1-pProcessBoth)*pTask1First proportion of trials. 
% - hSing: handle to the single-task points
% 
% by Alex White, 2017

function [pProcessBoth, pTask1First, hSing] = plotAOCWithPredictions(as,es,plotOpt)

if ~isfield(plotOpt,'markSz')
    plotOpt.markSz = 12;
end
if ~isfield(plotOpt,'doLegend')
    plotOpt.doLegend = true;
end
if ~isfield(plotOpt,'doAxLabels')
    plotOpt.doAxLabels = true;
end
if ~isfield(plotOpt,'plotSerialPrediction')
    plotOpt.plotSerialPrediction = true;
end

if ~isfield(plotOpt,'sideLabels')
    plotOpt.sideLabels = {'Left','Right'};
end

if ~isfield(plotOpt,'axLineWidth')
    plotOpt.axLineWidth = 1;
end

if ~isfield(plotOpt,'datLineWidth')
    plotOpt.datLineWidth = 1;
end

if ~isfield(plotOpt,'ebLineWidth')
    plotOpt.ebLineWidth = 1;
end

if ~isfield(plotOpt,'chanceLevel')
    plotOpt.chanceLevel = 0.5;
end

if ~isfield(plotOpt,'units')
    plotOpt.units = 'Ag';
end

if ~isfield(plotOpt,'edgeColors')
    plotOpt.edgeColors = repmat(hsv2rgb([0.6 0.7 0.4]), 2, 1);
end
if ~isfield(plotOpt,'fillColors')
    plotOpt.fillColors = [hsv2rgb([0.6 0.7 0.4]); ones(1,3)];
end

if ~isfield(plotOpt, 'predLineColor')
    plotOpt.predLineColor =  plotOpt.edgeColors(2,:);
end

%% setup axes
if ~isfield(plotOpt,'axlims')
    axlims = [plotOpt.chanceLevel 1];
else
    axlims = plotOpt.axlims;
end

if ~isfield(plotOpt,'axticks')
    axticks =  plotOpt.chanceLevel:(0.125):1;
else
    axticks = plotOpt.axticks;
end

if ~isfield(plotOpt,'axtickLabels')
    axtickLabels = cell(1,length(axticks));
    for ai=1:length(axticks)
        v = axticks(ai);
        if mod(ai,2)==1  %if ai==1 || ai==length(axticks)
            if (round(v)-v)==0
                frmt = '%.1f';
            elseif (round(10*v)-10*v) == 0
                frmt = '%.1f';
            elseif (round(100*v)-100*v) == 0
                frmt = '%.2f';
            else
                frmt = '%.3f';
            end
            axtickLabels{ai} = sprintf(frmt,v);
        else
            axtickLabels{ai} = '';
        end
    end
else
    axtickLabels = plotOpt.axtickLabels;
end
markSz = plotOpt.markSz;
axLineWidth = plotOpt.axLineWidth;
datLineWidth = plotOpt.datLineWidth;
ebLineWidth  = plotOpt.ebLineWidth;

%% set colors
singEdgeColor = plotOpt.edgeColors(1,:);
singFillColor = plotOpt.fillColors(1,:);
dualEdgeColor = plotOpt.edgeColors(2,:);
dualFillColor = plotOpt.fillColors(2,:);


%% Predictions of fixed-capacity parallel model:
if plotOpt.chanceLevel == 0.5
    singleAs = as(1,:); %data as area under ROC for single task conditions left and right
    %order needs flip:
    singleAs = fliplr(singleAs);
    singleDs = AgToDprime(singleAs); %convert to d'
    
    pSamplesOnTask1 = 0:0.001:1;  %vector of proportion of fixed number of sensory "samples" devited to task 1
    
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
    fixedA1s = DPrimeToAg(fixedD1s);
    fixedA2s =  DPrimeToAg(fixedD2s);
end
%% Plot the box for unlimited capacity model
hold on;

plot([as(1,2) as(1,2)],[axlims(1) as(1,1)],'k--','LineWidth',datLineWidth,'Color',plotOpt.predLineColor);
plot([axlims(1) as(1,2)], [as(1,1) as(1,1)],'k--','LineWidth',datLineWidth, 'Color',plotOpt.predLineColor);

%% Plot the serial model prediction of straight line:
plot([axlims(1) as(1,2)], [as(1,1) axlims(1)],'k-','LineWidth',datLineWidth, 'Color',plotOpt.predLineColor);

%% Plot fixed-capacity parallel prediction curve:
if plotOpt.chanceLevel == 0.5
    plot(fixedA1s,fixedA2s,'-','Color','k','LineWidth',datLineWidth, 'Color',plotOpt.predLineColor); %dualEdgeColor*0.8);
end
%% Solve for best-fitting serial model
[pProcessBoth, pTask1First, slope, intercept] = AnalyticDualTaskSerialModel(as(1,:),as(2,:), false);
%then predict dual-task performance given our 2 parameters:
predDual = SerialDualTaskAccGivenSingleTask(as(1,:),pTask1First,pProcessBoth);

%and plot it
%plot the best-fitting serial model line
if plotOpt.plotSerialPrediction
    
    leftSing = as(1,1);
    dRI = (leftSing-intercept)/slope; %predicted dual-task right side performance with 100% attention on left
    dR = as(1,2); %actual dual-task right side performance
    xs = [dRI dR];
    ys = slope*xs+intercept;
    plot(xs,ys,'r-');
end
%% Plot data
%Top single-task performance (side==1) on y axis
%add error bar:
plot(axlims([1 1]),as(1,1)+[-1 1]*es(1,1),'-','Color',singEdgeColor,'LineWidth',ebLineWidth);
hSing=plot(axlims(1),as(1,1),'o','MarkerSize',markSz,'MarkerFaceColor',singFillColor,'MarkerEdgeColor',singEdgeColor);

%bottom single-task performance (side==2) on x axis
%add error bar:
plot(as(1,2)+[-1 1]*es(1,2),axlims([1 1]),'-','Color',singEdgeColor,'LineWidth',ebLineWidth);
plot(as(1,2),axlims(1),'o','MarkerSize',markSz,'MarkerFaceColor',singFillColor,'MarkerEdgeColor',singEdgeColor);

%dual-task performance
%with error bars:
dualx = as(2,2); dualy=as(2,1);
dualxE = es(2,2); dualyE = es(2,1);
plot(ones(1,2)*dualx,dualy+[-1 1]*dualyE,'-','Color',dualEdgeColor,'LineWidth',ebLineWidth);
plot(dualx+[-1 1]*dualxE,ones(1,2)*dualy,'-','Color',dualEdgeColor,'LineWidth',ebLineWidth);

hDual = plot(dualx,dualy,'o','MarkerSize',markSz,'MarkerFaceColor',dualFillColor,'MarkerEdgeColor',dualEdgeColor,'LineWidth',datLineWidth+1);

%add predicted dual-task performance with best-fitting general serial model
if plotOpt.plotSerialPrediction
    plot(predDual(2), predDual(1),'r.','MarkerSize',6);
end


xlim(axlims);
ylim(axlims);
set(gca,'XTick',axticks,'XTickLabel',axtickLabels);
set(gca,'YTick',axticks,'YTickLabel',axtickLabels);
set(gca,'LineWidth',axLineWidth);
set(gca,'FontName','Helvetica')

if plotOpt.doAxLabels
    xlabel(sprintf('%s %s',plotOpt.sideLabels{2},plotOpt.units));
    ylabel(sprintf('%s %s',plotOpt.sideLabels{1},plotOpt.units ));
end
axis square;

if plotOpt.doLegend
    legend([hSing hDual],{'Single-task','Dual-task'},'Location','NorthWest'); legend boxoff;
end


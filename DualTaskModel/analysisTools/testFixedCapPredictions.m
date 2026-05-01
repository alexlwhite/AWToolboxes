%% Making fixed capacity parallel model prediction from accuracy measured as Area under the ROC curve (Ag)

%Input: single-task accuracy levels (Ag)
singleTaskAgs = [0.75 0.8]; %the 1st one is plotted on the y-axis; the 2nd on the x-axis (don't ask me why)

% convert to d'
singleTaskDs = sqrt(2)*norminv(singleTaskAgs);

%For making the fixed-capacity parallel model's prediction, 
% we simulate varying the proportion of a fixed number of sensory "samples"
% devited to task 1 (e.g., the left word). 
% Here's a vector from 0 to 1 of those proportions of samples: 
pSamplesOnTask1 = 0:0.001:1; 

% Now we use those to calculate predicted d's for each task, for each proportion of samples devoted to task 1. 
% How to calculate d' for diffent numbers of "samples"? 
% Given sample sizes ss1 and ss2: ss1 = ss2/2...
% std1 = sqrt(2)*std2
% generalized:
%  if ss1 = ss2*q, then
%  std1 = sqrt(1/q)*std2

%therefore, d' relation is:
% d1 = d2/sqrt(1/q)

%predicted dual-task dprimes
fixedD1s = singleTaskDs(1)./sqrt(1./pSamplesOnTask1);
fixedD2s = singleTaskDs(2)./sqrt(1./(1-pSamplesOnTask1));

%predicted dual-task Ags, converted from predicted dprime
predictedAgTask1 = normcdf(fixedD1s/sqrt(2));
predictedAgTask2 = normcdf(fixedD2s/sqrt(2));

%% plot the data points and the unlimited-capacity parallel  prediction box 
axlims = [0.5 1]; 
datLineWidth = 2;
plotOpt.predLineColor = 'k';

figure; hold on;

plot([singleTaskAgs(2) singleTaskAgs(2)],[axlims(1) singleTaskAgs(1)],'k--','LineWidth',datLineWidth,'Color',plotOpt.predLineColor);
plot([axlims(1) singleTaskAgs(2)], [singleTaskAgs(1) singleTaskAgs(1)],'k--','LineWidth',datLineWidth, 'Color',plotOpt.predLineColor);

%% Plot the serial model prediction of straight line:
plot([axlims(1) singleTaskAgs(2)], [singleTaskAgs(1) axlims(1)],'k-','LineWidth',datLineWidth, 'Color',plotOpt.predLineColor);

%% Plot fixed-capacity parallel prediction curve:
hAg = plot(predictedAgTask2,predictedAgTask1,'-','Color','k','LineWidth',datLineWidth, 'Color',plotOpt.predLineColor); %dualEdgeColor*0.8);


xlim(axlims); ylim(axlims);
axis square;


%% NOW let's do the same thing but assuming that our units are in proportion correct 

%Input: single-task accuracy levels (p(corr))
singleTaskPCs = [0.75 0.8]; %the 1st one is plotted on the y-axis; the 2nd on the x-axis (don't ask me why)

% convert to d': let's assume that proportion correct is for a one-sample,
% yes/no task. Then I believe the relationship is this: 
singleTaskDs = 2*norminv(singleTaskPCs);
%as opposed to sqrt(2)*norminv(singleTaskAgs)


%For making the fixed-capacity parallel model's prediction, 
% we simulate varying the proportion of a fixed number of sensory "samples"
% devited to task 1 (e.g., the left word). 
% Here's a vector from 0 to 1 of those proportions of samples: 
pSamplesOnTask1 = 0:0.001:1; 

%predicted dual-task dprimes
fixedD1s = singleTaskDs(1)./sqrt(1./pSamplesOnTask1);
fixedD2s = singleTaskDs(2)./sqrt(1./(1-pSamplesOnTask1));

%predicted dual-task Ags, converted from predicted dprime
%For a one-interval, yes/no proportion correct, the equation to convert
%from dprime is: 
% % PC = normcdf(dprime/2);
predictedPCTask1 = normcdf(fixedD1s/2);
predictedPCTask2 = normcdf(fixedD2s/2);

%% plot the data points and the unlimited-capacity parallel prediction box 
axlims = [0.5 1];
datLineWidth = 1;
plotOpt.predLineColor = 'r';

 hold on;

plot([singleTaskPCs(2) singleTaskPCs(2)],[axlims(1) singleTaskPCs(1)],'k--','LineWidth',datLineWidth,'Color',plotOpt.predLineColor);
plot([axlims(1) singleTaskPCs(2)], [singleTaskPCs(1) singleTaskPCs(1)],'k--','LineWidth',datLineWidth, 'Color',plotOpt.predLineColor);

%% Plot the serial model prediction of straight line:
plot([axlims(1) singleTaskPCs(2)], [singleTaskPCs(1) axlims(1)],'k-','LineWidth',datLineWidth, 'Color',plotOpt.predLineColor);

%% Plot fixed-capacity parallel prediction curve:
hPC = plot(predictedPCTask2,predictedPCTask1,'-','Color','k','LineWidth',datLineWidth, 'Color',plotOpt.predLineColor); %dualEdgeColor*0.8);


xlim(axlims); ylim(axlims);
axis square;

legh=legend([hAg hPC],{'Ag','p(corr)'})
title(legh, 'Prediction based on')
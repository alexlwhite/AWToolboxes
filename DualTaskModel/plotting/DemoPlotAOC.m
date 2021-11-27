%% Demo of plotting an AOC 

%% fake data: 

%mean proportion correct: 
%[single-task_1 single-task_2; 
% dual-task_1 dual-task_2]

pcs = [0.80 0.85; ...
       0.66 0.70]; 
  
%standard errors (for error bars):
sems = [0.02 0.03; ...
       0.05 0.03]; 
   
   
%% plot options 
opt.markSz = 10; 
opt.doLegend = true;
opt.doAxLabels = true;
opt.sideLabels = {'Left','Right'};
opt.axLineWidth = 1; 
opt.datLineWidth = 1; 
opt.units = 'p(correct)';
opt.edgeColors = repmat([0 0.7 0.3],2,1);
opt.fillColors = opt.edgeColors; 
opt.fillColors(2,:) = 1; 
opt.axticks = 0.5:0.1:1;

%whether to plot the prediction of the 'fit' generalized serial model: 
opt.plotSerialPrediction  = false;

%% plot
figure; 
[pProcessBoth, pTask1First] = plotAOCWithPredictions(pcs,sems,opt); 

fprintf(1,'\n\nThe generalized serial model estimates that both stimuli were processed on %.1f%% of dual-task trials\n', 100*pProcessBoth);
fprintf(1,'and that stimulus 1 was processed ''first'' on %.1f%% of dual-task trials\n', 100*pTask1First);

   
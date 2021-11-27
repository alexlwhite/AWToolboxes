%script to test dual-task model 
close all; clear;

%%Serial model

modelType = 1; %1=serial, 2=parallel

%Single-Task Ag values from the real data
singleTaskLeftAg = 0.82;  
singleTaskRightAg = 0.82;  


%Extra capacity parameter: on what proportion of dual-task trials both
%sides can be processed
pDualBoth = 0; 

%Congruency effects implemented by pFlip2Side2:
%On some proportion of trials, attend to side 2 only, but without
%realizing it, so report that side even if asked about the right 
%pFlip2Side2 ranges from -.5 (flip to side 1 half the time) to 0.5 (flip to side 2)
pFlip2Side2 = 0.2; 

%Bias towards one side in dual-task (on trials when only 1 is attended) 
%implemented by dualTaskAttnBias: 
%Ranges from 0 to 1, with 0 meaning never attend to the left, 0.5 being
%perfectly balanced, and 1 meaning always attend to the left
%We could assume that it's related to pFlip2Side2
% dualTaskAttnBias = 0.5-pFlip2Side2*2;
% dualTaskAttnBias(dualTaskAttnBias>1) = 1; 
% dualTaskAttnBias(dualTaskAttnBias<0) = 0; 

dualTaskAttnBias = 0.5;

if pFlip2Side2>0
    trueLeftAg = (singleTaskLeftAg - pFlip2Side2*0.5) / (1-pFlip2Side2);
    trueRightAg = singleTaskRightAg;
else
   trueRightAg = (singleTaskRightAg + pFlip2Side2*0.5) / (1+pFlip2Side2);
   trueLeftAg = singleTaskLeftAg;
end

attndLeftMean = AgToDprime(trueLeftAg); 
attndRightMean = AgToDprime(trueRightAg);

nSimTrials = 50000;

[Ag, AgValsByIndex, dualTaskAccCorr, corrValsByIndex, contingentAgs, contingentAgValsByIndex] = GeneralDualTaskModel(modelType,attndLeftMean, attndRightMean, pFlip2Side2, dualTaskAttnBias, pDualBoth, nSimTrials);

congruencyEffects = Ag(:,:,3)-Ag(:,:,2)

%Note: if pFlip2Side2 == 0, only negative correlations with 2 targets, at
%0 correlation with 0 targets. But if pFlip2Side2>0, negative correlation
%with >0 targets (especially for side2only) and positive corr with 0
%targets. Same for pFlip2Side2<0, but most negative corr for side1only.
%% Parallel model 
modelType = 2; %1=serial, 2=parallel

extraCapacity = 0; 

Ag  = GeneralDualTaskModel(modelType,attndLeftMean, attndRightMean, pFlip2Side2, dualTaskAttnBias, extraCapacity, nSimTrials);

congruencyEffects = Ag(:,:,3)-Ag(:,:,2)

%Note: the parallel model can also hae a complex pattern of negative and
%positive congruency effects if pFlip2Side2~=0.
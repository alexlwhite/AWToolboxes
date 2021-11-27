%% function [fixedPredAg1s, fixedPredAg2s] = computeFixedCapacityPrediction(singleTaskAgs)
% Compute the fixed-capacity parallel model's predction for dual-task
% accuracy in an Attention Operating Characteristic, given single-task
% accuracy levels 
%
% Inputs: 
% - singleTaskAgs: a 1x2 vector of single-task accuracy levels (Ag = area under roc),
% with the 1st element being the value that should be plotted along the X axis, and
% the 2nd being the value that should be plotted along the Y Axis 
% 
% Outputs: 
% - fixedPredAgs: 2x1000 matrix of predicted dual-task accuracy levels for
% 1000 different possible divisions of fixed processing resources between
% task 1 and task 2, ranging from 0% to 100% for task 1. Row 1 = Accuracy levels for 
% task 1 (corresponding to what is plotted on x axis); Row 2 = accuracy
% levels for task 2. 
% - pSamplesOnTask1: 1 1x1000 vector of proportion of samples on
% task/stimulus 1, ranging from 0 to 1. 

function [fixedPredAgs, pSamplesOnTask1] = computeFixedCapacityPrediction(singleTaskAgs)

singleTaskDs = AgToDprime(singleTaskAgs); %convert to d'


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
fixedPredD1s = singleTaskDs(1)./sqrt(1./pSamplesOnTask1);
fixedPredD2s = singleTaskDs(2)./sqrt(1./(1-pSamplesOnTask1));

%dual-task AROCs,
fixedPredAg1s = DPrimeToAg(fixedPredD1s);
fixedPredAg2s =  DPrimeToAg(fixedPredD2s);

fixedPredAgs = [fixedPredAg1s; fixedPredAg2s];
%% function [minDistToFixed, nearestFixedPSampTask1] = compareToFixedCapacityModel(as)
% For a given pattern of dual- and single-task accuracy levels, this
% function computes how far the dual-task point is from the nearest point
% on the prediction curve made by the fixed-capacity parallel model. Simply
% finds the Euclidean distance (minDistToFixed). It also finds the
% resource-sharing parameter (nearestFixedPSampTask1) that corresponds to that
% point. 
%
% Alex White, 2017
% 
%Inputs: 
%- as: 2x2 matrix of area-under-ROC-curve measures.
%    rows = single-task, dual-task
%    columns = left side, right side

%Outputs: 
% - minDistToFixed: distance between the dual-task data point and the nearest point on the 
%   fixed-capacity parallel curve. Negative if dual-task is worse than
%   predicted by the fixed-capacity parallel model (i.e. under the curve)
% - nearestFixedPSampTask1: for the nearest poin on the fixed-capacity
%   parallel curve, what is the proportion of "samples" devoted to task 1. 


function [minDistToFixed, nearestFixedPSampTask1] = compareToFixedCapacityModel(as)

%can't accept 100% correct 
as(as==1) = 0.999;

%% Predictions of fixed-capacity parallel model:
dualx = as(2,2); dualy=as(2,1);

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

%predicted dual-task AROCs,
fixedA1s = DPrimeToAg(fixedD1s);
fixedA2s =  DPrimeToAg(fixedD2s);

%% Find nearest fixed-capacity parallel point to data: 
fixedDists = sqrt((dualx-fixedA1s).^2 + (dualy-fixedA2s).^2);
minDistToFixed = min(fixedDists); 
pointI = find(fixedDists==minDistToFixed);
pointI = pointI(1);
nearestFixedPSampTask1 = pSamplesOnTask1(pointI);

%% set distance to be negative if dual-task performance is WORSE than predicted by fixed-capacity model 
nearestX = fixedA1s(pointI); nearestY = fixedA2s(pointI); 
%see which is farther from the origin 
nearestRad = sqrt(nearestX^2 + nearestY^2);

dataRad = sqrt(dualx^2 + dualy^2); 
if dataRad<nearestRad
    minDistToFixed = -1*minDistToFixed;
end

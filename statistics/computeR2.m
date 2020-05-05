%% function rSqr = computeR2(fitFun, xs, ys, fitParams);
% compute r^2 as the proportion of variance explained by a model fit to
% some data. 
% Inputs: 
% - fitFun: handle to the function that should be predicting the data
% - xs: x-values of the data 
% - ys: y-values of the data
% - fitParams: best fit parameters for fitFun, as a vector 
function rSqr = computeR2(fitFun, xs, ys, fitParams)

%y-values predicted by the function with these parameters: 
predYs  = fitFun(fitParams, xs);

%total sum of squares in the data
SSTot   = sum((ys - mean(ys)).^2);

%sum of squared residuals in the predicted values:
SSRes   = sum((ys - predYs).^2);

%r^2: proportion of variance explained 
rSqr    = 1 - SSRes/SSTot;


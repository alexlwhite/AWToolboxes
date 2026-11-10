%% function [betas, rSqr, r2adj, SSres, P] = linearRegressionWithStats(X, Y) 
% by Alex White, 2018
% Solve the linear model: Y = X*betas, given observered X and Y. 
% Returns the least-squares solution for betas, and also r-squared,
% adjusted r-squared, the sum of squared errors of the model fit, and the
% number of free parametres. 
% Inputs: 
% - X: a n x c matrix of predictor values 
% - Y: a n x j matrix of observed values 
% 
% Outputs: 
% - B: a c x j matrix of beta weights 
% - rSqr: the r-squared value for the fit: 1-SSres/SStot
% - rSqrAdj: r-squared adjusted by the degrees of freedom: 
%   1 - (SSres/(n-P)) / (SStot/(n-1))
%   where n is the total number of data points (numel(Y)) and P is the 
%   number of free parameters (numel(B), INCLUDING any intercept). 
%   For a single column of Y, this equals the familiar formula
%   1-(1-rSqr)*(n-1)/(n-p-1), where p is the number of predictors NOT 
%   counting the intercept. 
%   Note: rSqr and rSqrAdj are computed relative to the mean of Y, so they 
%   are only meaningful if X includes an intercept (a column of ones). 
%   Without an intercept, they can even be negative.
% - SSres: sum of squared errors between predicted values of Y and actual Y
% - P: number of free parameters (numel(B)), including any intercept. 

function [B, rSqr, rSqrAdj, SSres, P] = linearRegressionWithStats(X, Y) 

if size(X,1) ~= size(Y,1) 
    error('X and Y must have the same number of rows');
end

B = pinv(X)*Y;

%Compute r-squared
Yhat = X*B;
resids = Y - Yhat;
meanB = nanmean(Y(:));
SSres = sum(resids(:).^2);
SStot = sum((Y(:) - meanB).^2);
rSqr = 1-SSres/SStot;

%compute adjusted R2:
% rSqrAdj = 1 - (SSres/dfRes) / (SStot/dfTot)
% For a single column of Y and an X that includes an intercept column, this
% equals the familiar 1-(1-rSqr)*(N-1)/(N-p-1), where p is the number of
% predictors NOT counting the intercept.
nObs = numel(Y);   %total number of data points (all columns of Y)
P = numel(B);      %number of free parameters (all betas, including any intercept)
dfRes = nObs - P;  %residual degrees of freedom
dfTot = nObs - 1;  %SStot is computed around one grand mean

if dfRes>0 %only compute if there are more data points than free parameters
    rSqrAdj = 1 - (SSres/dfRes)/(SStot/dfTot);
else
    rSqrAdj = NaN;
end

function [slope, intercept, rSqr, corrR, corrP] = simpleLinearRegression(X, Y)

if ndims(X)>2 || (size(X,1)>1 && size(X,2)>1) 
    error('X must be a column vector'); 
end
if ndims(Y)>2 || (size(Y,1)>1 && size(Y,2)>1) 
    error('Y must be a column vector'); 
end

if size(X,1)==1 && size(X,2)>1
    X = X'; 
end
if size(Y,1)==1 && size(Y,2)>1
    Y = Y'; 
end

if length(X)~=length(Y)
    error('X and Y must be equal in length');
end
    
N = length(X);

%% Least-squares regression:
%D is a design matrix, with the first column being a constant offest
%(intercept), and the second column being the x-values that the slope will
%be multiplied with
D = [ones(N,1) X];

B = D\Y; %B is the beta weights, the solution to the equation P*B = Y; 

intercept = B(1);
slope = B(2);


%% compute R-squared
%y-values predicted by the model 
Y_hat = D*B;

resids =Y_hat - Y;
SSResid = sum(resids.^2);
SSTot = sum((Y - mean(Y)).^2);
rSqr = 1-SSResid/SSTot;

%% Correlation 
%Pearson: linear relationship between normally distributed continous variables 
[corrR, corrP] = corr(X,Y,'type','Pearson');
%Spearman: any monotomic relationship between continuous or ordinal variables, based on rank order
%[rhoS, pvalS] = corr(X,Y,'type','Spearman')


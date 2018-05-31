%% Generate fake data according to basic linear model: 
%Y = B0 + B1*X + e

N = 20;
X = randn(N,1); 
trueSlope = 0.5; 
trueInt = -.8; 

Y = trueSlope*X + trueInt + randn(N,1)*.3;

%% Least-squares regression:
P = [ones(N,1) X];

B = P\Y;

intercept = B(1);
slope = B(2);

%% Correlation 
%Pearson: linear relationship between normally distributed continous variables 
[rhoP, pvalP] = corr(X,Y,'type','Pearson') 
%Spearman: any monotomic relationship between continuous or ordinal variables, based on rank order
[rhoS, pvalS] = corr(X,Y,'type','Spearman')

%% compute R-squared
Y_hat = intercept + slope*X;

resids =Y_hat - Y;
SSResid = sum(resids.^2);
SSTot = sum((Y - mean(Y)).^2);
rSqr = 1-SSResid/SSTot;

%% Plot
figure; 
plot(X,Y,'b.'); 
hold on; 
plot([min(X) max(X)], [min(X) max(X)]*slope+intercept, 'k-'); 
xlabel('X');
ylabel('Y'); 

textX = min(X) + 0.02*(max(X)-min(X)); 
textY = min(Y) + 0.95*(max(Y)-min(Y)); 

text(textX, textY,sprintf('slope = %.3f, intercept = %.3f, r^2 = %.3f',slope, intercept, rSqr));

axis square; box off;
set(gcf,'color','w');


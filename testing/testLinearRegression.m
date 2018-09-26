%testing linear regression, in the context of 'inverted encoding model' or
%'forward model' 
% 
% Base model 
%y = aX1 + bX2 + error 
%for each of 3 conditions: 

a = [12 10 8];
b = [-0.5 -1 -2]; 
betas = [a; b]; 

n = 100; %number of data points
m = size(betas,1); %number of channels, one for X1 and one for X2
j = size(a,1); %number of conditions 

X = rand(n,m)*100; %predictors 

%true Y:
Y = X*betas;
%perturb by noise
Y = Y+randn(size(Y))*50;

[betasHat, rSqr, rSqrAdj, SSres, P] = linearRegressionWithStats(X, Y); 

%Try for each 'condition' separately 

betasHatSep = NaN(m,j);
rSqrAdjSep = NaN(1,j);
for ji=1:j
   y = X*betas(:,ji);
   y = y+randn(size(y))*50;

   [betasHatSep(:,ji), ~, rSqrAdjSep(ji)] = linearRegressionWithStats(X,y);
    
end
%turns out the same 
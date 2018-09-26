function [Fval, Pval] = nestedF(SSs, Ks, N)
%% function [Fval, Pval] = nestedF(SSs, Ks, N)
% by Alex White, 2018, using notes from Geoff Boynton 
%
%Computes a nested F test to compare two regresssion models, given sums of 
%squared errors for the two models, the numbers of parameters in the two models,
% and the total number of data points 
%
% Inputs: 
% - SSs: a 1x2 vector of sum of squared errors. SSs(1) = restricted
% (smaller) model; SSs(2) = complete (bigger) model. 
% - Ks: a 1x2 vector of number of parameters in each model, restricted and
% complete. Ks(1)<Ks(2)
% - N: number of data points 

%from Geoff's "Linear Regression and Anova" handout
%  F = ((SSr-SSc)/(kc-kr))/(SSc/(n-kc))
% where SSr is the sum of squared errors for the "restricted" model, the
% one with fewer regressors, SSc is the sum of squared errors for the
% 'complete' model, and kc is the total number of regressors in the
% complete model, kr is the number of regressors for the restricted model,
% and n is the total number of data points

SSr = SSs(1); 
SSc = SSs(2); 
kr = Ks(1); 
kc = Ks(2); 

if kr>kc
    error('(nestedF) Number of parameters for restrited model is bigger than for complete model!')
end

Fval = ((SSr-SSc)/(kc-kr))/(SSc/(N-kc));
df1 = kc-kr; 
df2 = N-kc;
Pval = 1-fcdf(Fval,df1,df2);

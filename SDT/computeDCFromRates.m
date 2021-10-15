function [d, c, c2, beta, hitR, FAR, corrected] = computeDCFromRates(hitR, FAR, nPres, nAbst, rateCorr)
%Compute signal detection theory variables from hit and false alarm rates and trial numbers
% by Alex White, 2014
% 
% Usage: [d, c, c2, beta, hitR, FAR, corrected] = computeDCFromRates(hitR, FAR, nPres, nAbs, rateCorr)
% Input: 
% - hitR is the hit rate, between 0 and 1
% - FAR is the false alarm rate, between 0 and 1
% - nPres is the number of trials with signal present 
% - nAbs is the number of trials with signal absent 
% - rateCorr is the amount by which to correct a hit or false alarm rate of
%   0 or 1
%
% Output: 
% - d, is dprime
% - c is criterion as distance from 0 (mean of noise distribution, depends
% only on false alarms
% - c2 is criterion as distance from the neutral point (d/2)
% - beta is criterion expresssed as the likelihood ratio of sensory evidence at c: (p(c | present) / p(c | absent)
% - hitR is the hit rate 
% - FAR is the false alarm rate 
% - corrected, whether or not a correction had to be applied to 100% hits
% or 0% false alarms. That's done by adding 1 error. 



corrected=false;

%Deal with case when there are too few trials, or hit rate is 1 or 0
if ~exist('rateCorr','var')
    rateCorrH = 1/(2*nPres);
    rateCorrF = 1/(2*nAbst);
else
    rateCorrH = rateCorr;
    rateCorrF = rateCorr;
end

if nPres<5
    if nPres>0
        if hitR==1, hitR=0.99;
        elseif hitR==0, hitR=0.01; end
        corrected=true;
    else
        hitR=NaN;
    end
   
elseif hitR==0
    hitR=rateCorrH;
    corrected=true;
elseif hitR==1
    hitR=1-rateCorrH;
    corrected=true;
end

if nAbst<5
    if nAbst>0
        if FAR==1, FAR=0.99;
        elseif FAR==0, FAR=0.01; end
        corrected=true;
    else
        FAR=NaN;
    end
elseif FAR==0
    FAR=rateCorrF;
    corrected=true;
elseif FAR==1
    FAR=1-rateCorrF;
    corrected=true;
end


d=norminv(hitR)-norminv(FAR);

%Criterion: 
%First, just the z-score of the correct rejection rate. Simple, on same
%units as d'
c=norminv(1-FAR); 

%Second formula: a measure of how far the criterion is from neutral, which
%is d/2. This is equivalent to c-d/2 (using the first c)
c2= -0.5*(norminv(hitR) + norminv(FAR));

%Third, beta: the likelihood ratio at c 
beta = normpdf(c,d,1)/normpdf(c,0,1); 

%other transformations relevant to these criterion measues: 
%if you expand normdf in the above formula for beta: 
%beta = exp(c*d - 0.5*d^2)
% 
% We can also say that: 
%c2 = c - d/2
% 
% So we can also express beta in terms of c2: 
% 
% beta  = exp(c2*d) 
% 
% That means that 
% c2 = log(beta)/d

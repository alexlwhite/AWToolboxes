function [d, c, c2, beta, hitR, FAR, corrected] = computeDC(pres, resp)
%Compute d' (d) and criterion (c) 
% by Alex White, 2011
% 
% Usage: [d, c, c2, corrected] = computeDC(pres, resp)
% Input: 
% - pres should be vector of 1s (for signal present) and 0s (for signal
% absent) 
% 
% - resp should be vector of 1s (for 'yes' response) and 0s (for 'no'
% response)
%
% Output: 
% - d, is dprime
% - c is criterion as distance from 0 (mean of noise distribution, depends
% only on false alarms
% - c2 is criterion as distance from the neutral point (d/2)
% - beta is the likelihood ratio at c
% - hitR is the hit rate 
% - FAR is the false alarm rate 
% - corrected, whether or not a correction had to be applied to 100% hits
% or 0% false alarms. That's done by adjusting the rate to what it would be if he 
% had run twice as many trials and got 1 different response. 


hits=(pres & resp); 
FAs=(~pres & resp); 

hitR=sum(hits)/sum(pres); 
FAR=sum(FAs)/sum(~pres); 

corrected=false;

%Deal with case when there are too few trials, or hit rate is 1 or 0
if sum(pres)<2
    if sum(pres)>0
        if hitR==1, hitR=0.99;
        elseif hitR==0, hitR=0.01; end
        corrected=true;
    else
        hitR=NaN;
    end
    
elseif hitR==0
    hitR=1/(2*sum(pres));
    corrected=true;
elseif hitR==1
    hitR=(2*sum(pres)-1)/(2*sum(pres));
    corrected=true;
end

%Deal with case when there are too few trials, or FA rate is 1 or 0
if sum(~pres)<2
    if sum(~pres)>0
        if FAR==1, FAR=0.99;
        elseif FAR==0, FAR=0.01; end
        corrected=true;
    else
        FAR=NaN;
    end
elseif FAR==0
    FAR=1/(2*sum(~pres));
    corrected=true;
elseif FAR==1
    FAR=(2*sum(~pres)-1)/(2*sum(~pres));
    corrected=true;
end


d=norminv(hitR)-norminv(FAR);

%don't accept d' less than 0. 
%if d<0, d=0; end;

%Criterion: 
%First, just the z-score of the correct rejection rate. Simple, on same
%units as d'
c=norminv(1-FAR); 

%Second formula: a measure of how far the criterion is from neutral, which
%is d/2. This is equivalent to c-d/2 (using the first c)
c2= -0.5*(norminv(hitR) + norminv(FAR));

%Third, beta: the likelihood ratio at c  
beta = normpdf(c,d,1)/normpdf(c,0,1); 


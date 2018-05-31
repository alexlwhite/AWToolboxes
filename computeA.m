%% function [A, b] = computeA(h,f) 
% Compute sensitivity meausure "A" and bias measure "b", as defined by
% Zhang & Mueller, Psychometrika, 2005
%
% Inputs
% - h: hit rate
% - f: false alarm rate 
% 
% Outputs: 
% - A: sensitivity measure. It is the average areas of minimum and maximum
% ROC curves that pass through the single point (f,h). Ranges from 0.5
% (chance) to 1 (perfect).
% - b: a measure of bias related to the slope of the ROC curve at that
% point 


function [A, b] = computeA(h,f) 

if (f<=.5 && h>=.5)
    A = .75 + (h-f)/4 - f*(1-h);
elseif (f<=h && h<.5)
    A = .75 + (h-f)/4 - f/(4*h);
else 
    A = .75 + (h-f)/4 - (1-h)/(4 * (1-f));
end


if (f<=.5 && h>=.5)
    b = (5-4*h)/(1+4*f);
elseif (f<h && h<.5)
    b = (h^2+h)/(h^2+f);
else 
    b =((1-f)^2 + (1-h))/((1-f)^2 + (1-f));
end

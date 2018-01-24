% function [HR, FR] = computeROCRates(pres, resp) 
% Alex White, April 2015
% 
% Computes hit rate (HR) and false alarm rate (FAR) at multiple decision
% "criteria" from detection task responses, wherein the observer has
% multiple keys to report target presence/absence at varying levels of
% "confidence." 
% 
% Inputs: 
% - pres: a vector of 1s and 0s, with 1 entry for each trial indicating
% whether the target stimulus was present (1) or absent (0).
% - resp: a vector that indicates the observer's response on each trial
% (corresponding to the trials in pres). The elements of resp should be
% numbers that indicate which key in a set of K keys was pressed, ranked 
% from low to high confidence of target presence (e.g., 1=sure no, 2=maybe no,
% 3=maybe yes, 4=sure yes). 
% - keys: a vector of identifiers for the possible keys 
% 
% This function computes hit rate and false alarm rate at K-1 criteria
% corresponding to positions between the possible response keys. These can
% then be used to compute the area under the ROC. 
% 
% Outputs: 
% - HR: a vector (length K-1), of hit rates at all possible
% criteria. The order is from conservative to liberal (e.g., the first
% element of HR includes only responses of the highest confidence, and the last
% element includes responses for all levels of confidence but the first
% (which at that criterion corresponds to report of target absence). 
% - FR: a vector (length K-1) of false alarm rates, organized just like HR.
% 

function [HR, FR] = computeROCRates(pres, resp, keys) 

%number of possible responses:
nKeys = length(keys); 

%another way: 
HR = zeros(1,nKeys);
FR = zeros(1,nKeys);
keyOrder = nKeys:(-1):1;
for i = 1:length(keyOrder)
    keyI = keyOrder(i);
    HR(i) = nanmean(resp(pres == 1)>=keys(keyI));
    FR(i) = nanmean(resp(pres == 0)>=keys(keyI));
end

%add the lower right corner, which would correspond to responses greater
%than the highest posstible response, but anchors the ROC curve 
HR = [0 HR];
FR = [0 FR];    

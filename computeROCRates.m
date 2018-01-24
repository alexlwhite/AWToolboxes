% function [HR, FR] = computeROCRates(pres, resp, keys) 
% Alex White, April 2015
% 
% Computes hit rate (HR) and false alarm rate (FAR) at multiple decision
% "criteria" from detection task responses, wherein the observer has
% multiple keys to report target presence/absence at varying levels of
% "confidence." 
% 
% Inputs: 
% - pres: a vector of 1s and 0s, with 1 entry for each trial indicating
%   whether the target stimulus was present (1) or absent (0).
% - resp: a vector that indicates the observer's response on each trial
%   (corresponding to the trials in pres). The elements of resp should be
%   numbers that indicate which key in a set of K keys was pressed, ranked 
%   from low to high confidence of target presence (e.g., 1=sure no, 2=maybe no,
%    3=maybe yes, 4=sure yes). 
% - keys: a vector of identifiers for the possible keys. If this is left
%   out, keys = unique(resp)
% 
% This function computes hit rate and false alarm rate at K-1 criteria
% corresponding to positions between the possible response keys. These can
% then be used to compute the area under the ROC. 
% 
% Outputs: 
% - HR: a vector (length K-1), of hit rates at all possible
%   criteria. The order is from conservative to liberal (e.g., the first
%   element of HR includes only responses of the highest confidence, and the last
%   element includes responses for all levels of confidence but the first
%   (which at that criterion corresponds to report of target absence). 
% - FR: a vector (length K-1) of false alarm rates, organized just like HR.
% 

function [HR, FR] = computeROCRates(pres, resp, keys) 

if nargin<3
    %response key identifiers: 
    keys = unique(resp);
end

%number of possible responses:
nKeys = length(keys); 

%proportions of target-present trials in which a particular key was pressed
P = zeros(1,nKeys);
%proportions of target-absent trials in which a particular key was pressed
A = zeros(1,nKeys);

presTrials = pres == 1;
abstTrials = pres == 0;

for keyI = 1:nKeys
    P(keyI) = nanmean(resp(presTrials)==keys(keyI));
    A(keyI) = nanmean(resp(abstTrials)==keys(keyI));
end

%Compute HR and FR with cumulative sum of response probabilities P and A. 
% But reverse the order before doing the cumulative sum, because the elements on
%P and A are organized from lowest to highest confidence, and HR and FR go
%from conservative to liberal. 
HR = cumsum(P(nKeys:-1:2));
FR = cumsum(A(nKeys:-1:2));
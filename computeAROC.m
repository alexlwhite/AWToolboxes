% function A = computeAROC(HR, FR) 
% Alex White, April 2015
% 
% This function estimates the area under the ROC curve (A) for data from a
% detection task with reports of varying levels of confidence. 
% 
% Inputs: 
% - HR: a vector of length C, containing the hit rates compute at C
% different criteria, from conservative to liberal. HR should therefore be
% monotonically increasing. It's best if HR includes the most extreme
% criteria that produce hit rates of 0 and 1, to anchor the curve at the
% corners (even if that wasn't really possible for the subject to report in
% the experiment). The corresponding function computeROCRates should do
% that. But if the 1st element of HR is greater than 0, we'll concatenate it
% with a 0, and if the last element of HR is less than 1, we'll concatenate
% a 1 onto it. 
% - FR: a vector (also of length C) containing the corresponding false
% alarm rates. Similarly, it should begin at 0  and to go 1. If not, those
% value are concatenated. 
% 
% Outputs: 
% - A: the estimated area under the ROC curve made of points (FR, HR). 
% This is computed by summing the areas of rectangles formed below the points 
% and right triangles formed between them. The points (0,0) and (1,1) 
% are added to complete the curve. 


function A = computeAROC(HR, FR) 


%add the absolute corners if necessary:
if HR(1)>0
    HR = [0 HR]; 
end
if HR(end)<1
    HR = [HR 1];
end
if FR(1)>0
    FR = [0 FR]; 
end
if HR(end)<1
    FR = [FR 1];
end


%compute area of 'rectangles' 
nRect = length(FR)-2;
rectA = zeros(1,nRect);
for ri = 1:nRect
    rectA(ri) = (FR(ri+2)-FR(ri+1))*HR(ri+1); 
end

%compute area of 'trianges'
nTri = nRect+1;
triA = zeros(1,nTri); 
for ti = 1:nTri
   triA(ti) = 0.5*(FR(ti+1)-FR(ti))*(HR(ti+1)-HR(ti));  
end

A = sum(rectA)+sum(triA);


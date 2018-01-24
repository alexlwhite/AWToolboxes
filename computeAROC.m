% function A = computeAROC(HR, FR) 
% Alex White, April 2015
% 
% This function estimates the area under the ROC curve (A) for data from a
% detection task with reports of varying levels of confidence. 
% 
% Inputs: 
% - HR: a vector of length C, containing the hit rates compute at C
% different criteria, from conservative to liberal. HR should therefore be
% monotonically increasing. 
% - FR: a vector (also of length C) containing the corresponding false alarm rates. 
% 
% Outputs: 
% - A: the estimated area under the ROC curve made of points (FR, HR). 
% This is computed by summing the areas of rectangles formed below the points 
% and right triangles formed between them. The points (0,0) and (1,1) 
% are added to complete the curve. 


function A = computeAROC(HR, FR) 


%add the absolute corners: 
allh = [0 HR 1]; 
allf = [0 FR 1]; 


%compute area of 'rectangles' 
nRect = length(FR);
rectA = zeros(1,nRect);
for ri = 1:nRect
    rectA(ri) = (allf(ri+2)-allf(ri+1))*allh(ri+1); 
end

%compute area of 'trianges'
nTri = nRect+1;
triA = zeros(1,nTri); 
for ti = 1:nTri
   triA(ti) = 0.5*(allf(ti+1)-allf(ti))*(allh(ti+1)-allh(ti));  
end

A = sum(rectA)+sum(triA);


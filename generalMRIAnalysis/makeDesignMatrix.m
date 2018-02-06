%% function X = makeDesignMatrix(numCond, numTR, lengthHRF, eventSequence, doPlot)
% Make a design matrix for event-related fMRI analysis 
% 
% Inputs: 
% - eventSequence: a Tx1 vector of integers that describes which conditions
%   were on at which time points. T is the total number of TRs in the scan
% - conds: a vector of integers specifying which event types to pull out of eventSequence. 
%   Simple case: if you have 3 conditions in eventSequence, then conds = 1:3. 
% - lengthHRF: the length of the HRF that will be convolved with this
%   design matrix to make predicted responses
% - doPlot: whether to show a grascale image of the design matrix 
% 
% Outputs: 
% - X: the design matrix, of 0s and 1s. 
%   number of rows = numTR. 
%   number of columns = number of conditions x lengthHRF 
%   The first column has 1s for each time point when the 1st
%   condition was on. The 2nd column is the same but shifted down by 1 TR,
%   and so on, such that the 1st lengthHRF columns represent condition 1,
%   and then the next lengthHRF columns represent condition 2, etc. 


function X = makeDesignMatrix(eventSequence, conds, lengthHRF, doPlot)

X = []; 

numTR = length(eventSequence);

for j=conds
    Xj = zeros(numTR,lengthHRF);
    temp = eventSequence == j;
    
    for i=1:lengthHRF
        Xj(:,i) = temp;
        temp = [0;temp(1:end-1)];
    end
    X = [X,Xj];
end

if doPlot
    figure;  
    imagesc(X)
    axis equal
    colormap(gray)
    axis off
    title('Design matrix');
end
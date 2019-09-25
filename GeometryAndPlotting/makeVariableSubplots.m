%% function subplotPositions = makeVariableSubplots(nRows, nCols, margins, verticalSpace, horizontalSpace, relativeHeights, relativeWidths)
% Create a matrix of subplot coordinates for a figure with multiple panels.
% The coordinates can be used with subplot('postion', [left bottom width height])
%  
% Inputs: 
% - nRows: number of rows in the figure; 
% - nCols: number of columns in the figure; 
% - margins: a 1x4 vector of margin space [top bottom left right],
%   expressed as proportions of the whole figure dimensions. 
% - verticalSpace: how much space between rows. Either a single number, or
%   if that spacing should vary between rows, then this variable should be a vector of size 1x(nRows-1)  
% - horizontalSpace: similar to verticalSpace, but for the horizontal space
%   between columns. Either a single number, or a vector of size 1x(nCols-1)
% - relativeHeights: relative heights of the rows. Of size 1xnRows. For
%   instance, if there are 3 rows, to make the top row be twice as tall as
%   the other two, relativeHeights = [2 1 1]. 
% - relativeWidths: similar to relativeHeights, except for each column's
%   width. 
%
% Output:   
% - subplotPositions: a matrix of coordinates of size nRows x nCols x 4. 
%   For the subplot in the 2nd row and 3rd column,  then you would pull out
%   subplotPositions(2,3,:). Those 4 numbers are [leftEdge bottomEdge width height] of the subplot, 
%   to be used in the subplot('position') function. 
% 
% by Alex L. White, at the University of Washington, 2019 
% 
function subplotPositions = makeVariableSubplots(nRows, nCols, margins, verticalSpace, horizontalSpace, relativeHeights, relativeWidths)

topMargin    = margins(1); 
bottomMargin = margins(2); 
leftMargin   = margins(3); 
rightMargin  = margins(4); 


if length(verticalSpace)==1
    verticalSpace = ones(1,nRows-1)*verticalSpace;
end
if length(horizontalSpace)==1
    horizontalSpace = ones(1,nCols-1)*horizontalSpace;
end

availableHeight =  (1-topMargin-bottomMargin-sum(verticalSpace));
subplotHeights = availableHeight*relativeHeights/sum(relativeHeights);

availableWidth = (1-leftMargin-rightMargin-sum(horizontalSpace));
subplotWidths = availableWidth*relativeWidths/sum(relativeWidths);

subplotPositions = zeros(nRows,nCols,4);

for ri=1:nRows
    for ci=1:nCols
        leftPos = leftMargin + sum(horizontalSpace(1:(ci-1))) + sum(subplotWidths(1:(ci-1)));
        
        plotsBelow = (ri+1):nRows;

        bottomPos = bottomMargin +  sum(verticalSpace(plotsBelow-1)) + sum(subplotHeights(plotsBelow));
        
        subplotPositions(ri,ci,:) = [leftPos bottomPos subplotWidths(ci) subplotHeights(ri)];
    end
end
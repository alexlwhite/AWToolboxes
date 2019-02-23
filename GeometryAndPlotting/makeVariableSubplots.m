function subplotPositions = makeVariableSubplots(nRows, nCols, margins, verticalSpace, horizontalSpace, relativeHeights, relativeWidths)

topMargin    = margins(1); 
bottomMargin = margins(2); 
leftMargin   = margins(3); 
rightMargin  = margins(4); 


if length(verticalSpace)==1
    verticalSpace = ones(1,nRows-1)*verticalSpace;
end
if length(horizontalSpace)==1
    horizontalSpace = ones(nCols-1)*horizontalSpace;
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
nRows = 3;
nCols = 3; 

nC = nRows*nCols;

figure;
for c = 1:nC
    subplot(3,3,c); 
    text(1,1,sprintf('%i',c)); 
    
    [colI, rowI] = ind2sub([nCols nRows], c);
    
    text(0, 0, sprintf('row=%i, col=%i', rowI, colI));
    
    xlim([-2 2]); ylim([-2 2]);
end
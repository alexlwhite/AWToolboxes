nRows = 8;
nCols = 8; 

nC = nRows*nCols;

figure;
for c = 1:nC
    
    subplot(nRows, nCols, c); 
    text(0,0,sprintf('%i',c)); 
    
    [colI, rowI] = ind2sub([nCols nRows], c);
    
    %text(0, 0, sprintf('row=%i, col=%i', rowI, colI));
    
    xlim([-2 2]); ylim([-2 2]);
    axis off;
end

%% 
figure; 
subplot(nRows, nCols, [1 5 9]);
text(0, 0, '[1 5 9]');


figure;
nrows = 4;
ncols = 2;
for rowi=1:nrows
    for coli=1:ncols
        subi=sub2ind([ncols nrows],coli,rowi); 
        subplot(nrows,ncols,subi);
        text(0,0,sprintf('row=%i,col=%i,subi=%i',rowi,coli,subi), 'HorizontalAlignment','center');
        xlim([-3 3]); ylim([-2 2]);
    end
end
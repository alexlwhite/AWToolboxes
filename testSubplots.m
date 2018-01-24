figure;
nrows = 4;
ncols = 2;
for rowi=1:nrows
    for coli=1:ncols
        subi=sub2ind([ncols nrows],coli,rowi); 
        subplot(nrows,ncols,subi);
        text(1,1,sprintf('row=%i,col=%i,subi=%i',rowi,coli,subi));
        xlim([-3 3]); ylim([-2 2]);
    end
end
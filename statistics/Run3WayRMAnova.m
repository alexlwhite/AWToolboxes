%Run3WayRMAnova(ds,file)
% Runs a 3-way repeated measures ANOVA on input matrix ds, using the
% RMAOV33Save function. 
% 
%Inputs: 
% - ds, a 4-dimensional matrix first 3 dimensions for each of 3 factors, 4th dimension for subjects 
% - file: handle to a file where the results table should be printed. Input
%   1 to print to command window. 

function Run3WayRMAnova(ds,file)

nsubj = size(ds,ndims(ds));

ndat = numel(ds); 
y = NaN(1,ndat); 
g1 = y; 
g2 = y;
g3 = y; 
g4 = y; 

dpnum=0; gnum=0;
for v1=1:size(ds,1)
    
    gnum=gnum+1;
    for v2=1:size(ds,2)
        for v3=1:size(ds,3)
            for vsubj=1:nsubj
                dpnum=dpnum+1;
                y(dpnum)=ds(v1,v2,v3,vsubj);
                g1(dpnum)=v1;
                g2(dpnum)=v2;
                g3(dpnum)=v3;
                g4(dpnum)=vsubj;
            end
        end
    end
end



dataMat = [y' g1' g2' g3' g4'];

RMAOV33Save(dataMat,0.05,file);
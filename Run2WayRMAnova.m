%stats = Run2WayRMAnova(ds, file) 
%
% Runs a repeated-measures two-way anova on data in matrix "ds"
% uses the function rm_anova2, taken from the Matlab exchange 
% 
% Inputs: 
% - ds: matrix of data values with 3 dimensions. 1st dimension is
% experimental factor 1; 2nd dimension is experimental factor 2; 3rd
% dimension is subjects 
% - file: handle to a file to print the ANOVA table to. Input "1" to print
% to command line 
% - factnames: a 1x2 cell array with names of the 2 factors 
% 
% Outputs: 
% - stats: a cell array containing the full ANOVA statistical table. This
% the same as what will be printed to "file". 

function stats = Run2WayRMAnova(ds, file, factnames)

g1=zeros(1,numel(ds));
g2=zeros(1,numel(ds));
g3=zeros(1,numel(ds));
y=zeros(1,numel(ds));

dpnum=0;
for v1=1:size(ds,1)
    for v2=1:size(ds,2)
        for v3=1:size(ds,3)
            dpnum=dpnum+1;
            g1(dpnum)=v1;
            g2(dpnum)=v2;
            g3(dpnum)=v3;
            
            y(dpnum)=ds(v1,v2,v3);
        end
    end
end

stats=rm_anova2(y,g3,g1,g2,factnames);

for sxi=1:size(stats,1)
    for syi=1:size(stats,2)
        if ischar(stats{sxi,syi})
            fprintf(file,sprintf('%s\t',stats{sxi,syi}));
        else
            fprintf(file,sprintf('%.6f\t',stats{sxi,syi}));
        end
    end
    fprintf(file,sprintf('\n'));
end
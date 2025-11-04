function compRes = printLMEModelComparison(smallModel, complexModel, modelNames, sf, checkNesting)

if nargin<5
    checkNesting = true;
end

%model comparison
compRes = compare(smallModel, complexModel,'CheckNesting',checkNesting);

%print it out
fprintf(sf, '\nComparing two LME model fits: %s vs %s\n', modelNames{1}, modelNames{2});
cProps = compRes.Properties.VarNames;
fprintf(sf,'\tModelName');
for cpi = 1:length(cProps)
    fprintf(sf, '\t%s', cProps{cpi});
end
fprintf(sf,'\n\t');
for mdi = 1:2
    fprintf(sf,'%s\t', modelNames{mdi});
    for cpi = 1:length(cProps)
        eval(sprintf('dat = compRes.%s(%i);', cProps{cpi}, mdi));
        fprintf(sf, '%.6f\t', dat);
    end
    fprintf(sf,'\n\t');

end
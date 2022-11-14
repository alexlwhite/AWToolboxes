function printLmeRes(lme, sf)

cs = lme.Coefficients;
fprintf(sf,'LME: %s', lme.Formula);
fprintf(sf,'\nFactor\tMean estimate\tt\tDF\tp');
for ci=1:length(cs.Name)
    fprintf(sf,'\n%s\t%.4f\t%.4f\t%i\t', cs.Name{ci}, cs.Estimate(ci), cs.tStat(ci), cs.DF(ci));
    if cs.pValue(ci)<0.001
         fprintf(sf,'%.2e', cs.pValue(ci));
    else
        fprintf(sf, '%.5f', cs.pValue(ci));
    end

end

anv = lme.anova;
fprintf(sf,'\n\nANOVA');
fprintf(sf,'\nFactor\tF\tDF1\tDF2\tp');
for ci=1:length(anv.Term)
    fprintf(sf,'\n%s\t%.4f\t%i\t%i\t', anv.Term{ci}, anv.FStat(ci), anv.DF1(ci), anv.DF2(ci));
    if anv.pValue(ci)<0.001
         fprintf(sf,'%.2e', anv.pValue(ci));
    else
        fprintf(sf, '%.5f', anv.pValue(ci));
    end
end
function printLmeRes(lme, sf)

cs = lme.Coefficients;
fprintf(sf,'LME: %s', lme.Formula);
fprintf(sf,'\nFactor\tMean estimate\tt\tDF\tp');
for ci=1:length(cs.Name)
    fprintf(sf,'\n%s\t%.4f\t%.4f\t%.4f\t%.7f\t', cs.Name{ci}, cs.Estimate(ci), cs.tStat(ci), cs.DF(ci), cs.pValue(ci));
end

anv = lme.anova;
fprintf(sf,'\n\nANOVA');
fprintf(sf,'\nFactor\tF\tDF1\tDF2\tp');
for ci=1:length(anv.Term)
    fprintf(sf,'\n%s\t%.4f\t%.4f\t%.4f\t%.7f\t', anv.Term{ci}, anv.FStat(ci), anv.DF1(ci), anv.DF2(ci), anv.pValue(ci));
end
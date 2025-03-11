%% function printLmeRes(lme, sf, printANOVA)
% prints to a text file the results of a linear mixed effects model fit. 
% Inputs: 
% - lme: the object output by fitlme or fitglme
% - sf: a handle to the file to print to (if 1, then it prints to command
% window) 
% - printANOVA: true or false, whether to also print the ANOVA-formatted
% results. 
%
% By Alex White, at Barnard College 

function printLmeRes(lme, sf, printANOVA)

if nargin<3
    printANOVA = true;
end

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
fprintf(sf, '\n');


if printANOVA
    anv = lme.anova;
    fprintf(sf,'\nANOVA');
    fprintf(sf,'\nFactor\tF\tDF1\tDF2\tp');
    for ci=1:length(anv.Term)
        fprintf(sf,'\n%s\t%.4f\t%i\t%i\t', anv.Term{ci}, anv.FStat(ci), anv.DF1(ci), anv.DF2(ci));
        if anv.pValue(ci)<0.001
            fprintf(sf,'%.2e', anv.pValue(ci));
        else
            fprintf(sf, '%.5f', anv.pValue(ci));
        end
    end
end
fprintf(sf, '\n');

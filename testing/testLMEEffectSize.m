meanGroupDiff = 0.15; 
meanCondDiff = 0.25; 

grandMean = 0.6; 
subjSigma = 0.2; 

N = 100; 


G1C1 = randn(N,1)*subjSigma + grandMean - meanGroupDiff/2 - meanCondDiff/2;
G1C2 = randn(N,1)*subjSigma + grandMean - meanGroupDiff/2 + meanCondDiff/2;
G2C1 = randn(N,1)*subjSigma + grandMean + meanGroupDiff/2 - meanCondDiff/2;
G2C2 = randn(N,1)*subjSigma + grandMean + meanGroupDiff/2 + meanCondDiff/2;

avgs = [mean(G1C1) mean(G1C2); mean(G2C1) mean(G2C2)]


allYs = {G1C1, G1C2; G2C1, G2C2};

D = table;

nRows = N*2*2; 

groupLabs = {'Group1','Group2'};
condLabs = {'Cond1','Cond2'};
groups = cell(nRows,1); 
conds = cell(nRows,1); 
subjects = NaN(nRows,1);
ys = NaN(nRows,1);

for gi=1:2
    for ci=1:2
       rowIs = (1:N) + (gi-1)*2*N + (ci-1)*N;
       subjects(rowIs) = (1:N) + (gi-1)*N;
       groups(rowIs) = groupLabs(gi);
       conds(rowIs) = condLabs(ci); 
       ys(rowIs) = allYs{gi,ci};
    end
end

D.y = ys; 
D.group = groups; 
D.cond = conds; 
D.subject = subjects;

eqtn = 'y~group*cond + (1|subject)'; 
lme = fitlme(D, eqtn, 'DummyVarCoding','effects'); 
display(lme); 
display(lme.anova);


lme = fitlme(D, eqtn, 'DummyVarCoding','reference'); 
display(lme); 
display(lme.anova);
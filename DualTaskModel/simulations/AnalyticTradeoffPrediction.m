%% try an analytic solution to predict stimulsu processing tradeoff according to the serial model 
% this assumes that when a side is not processed, p(correct) = 0.5. Simple
% guessing rule. 
% It's the pure all-or-none model. It never happens that both sides are
% processed. 
% Aug 9 2019: I don't quite understand why this predicts bigger effects
% than my simulations in which the unprocessed stimulus is represented by a
% draw from the 'default' distribution midway between target-absent and
% target-present. 

clear;

pStim1First = 0.5; %no bias to one side or the other

singleTaskDprimes = 0:0.2:4;
nDs = length(singleTaskDprimes);

pCorrect_OtherSideCorrect = NaN(nDs,2);
pCorrect_OtherSideIncorrect = NaN(nDs,2);

for di = 1:nDs
    attnd1Mean = singleTaskDprimes(di);
    attnd2Mean = attnd1Mean;
    
    singleAgs = DPrimeToAg([attnd1Mean attnd2Mean]);

    for si=1:2
        
        %other side processed
        if si==1
            pProcess = pStim1First;
        else
            pProcess = 1-pStim1First;
        end

        pProcessOther = 1-pProcess;
        
        %Ag | other side is 'processed' = 0.5;
        %Ag | this side is processed = singleAgs(si)
        
        %Ag | other side is correct   = [p(c) | (other side is processed & correct)] + [p(c) | (other side is not processed & correct)]
        %Ag | other side is incorrect = [p(c) | (other side is processed & incorrect)] + [p(c) | (other side is not processed & incorrect)]
        
        pOtherCorrect = pProcessOther*singleAgs(3-si) + 0.5*(1-pProcessOther);

        pOtherSideProcessedAndCorrect = singleAgs(3-si)*pProcessOther;
        pOtherSideProcessedAndIncorrect = (1-singleAgs(3-si))*pProcessOther;
        pOtherSideNotProcessedAndCorrect = 0.5*(1-pProcessOther);
        pOtherSideNotProcessedAndIncorrect = 0.5*(1-pProcessOther);
        
        pCorrect_OtherSideProcessedAndCorrect   = 0.5*pOtherSideProcessedAndCorrect;
        pCorrect_OtherSideProcessedAndIncorrect = 0.5*pOtherSideProcessedAndIncorrect;

        pCorrect_OtherSideNotProcessedAndCorrect   = singleAgs(si)*pOtherSideNotProcessedAndCorrect;
        pCorrect_OtherSideNotProcessedAndIncorrect = singleAgs(si)*pOtherSideNotProcessedAndIncorrect;
        

        pCorrect_OtherSideCorrect(di,si) = (pCorrect_OtherSideProcessedAndCorrect + pCorrect_OtherSideNotProcessedAndCorrect)/pOtherCorrect;
        pCorrect_OtherSideIncorrect(di,si) = (pCorrect_OtherSideProcessedAndIncorrect + pCorrect_OtherSideNotProcessedAndIncorrect)/(1-pOtherCorrect);
        
    end
end
    
%average over sides
tradeoffAccs = [mean(pCorrect_OtherSideIncorrect,2) mean(pCorrect_OtherSideCorrect,2)];

%% plot
hues = [0.33 0.67];
sats = [0.8 0.8];
vals = [0.7 1];

colrs = hsv2rgb([hues' sats' vals']);

figure; hold on;
plot([0.5 1], [0.5 1],'k--');
plot(tradeoffAccs(:,1), tradeoffAccs(:,2), '.-', 'Color',colrs(1,:));


xlabel('Ag | other side incorrect');
ylabel('Ag | other side correct');

%xlim([0.5 1]); ylim([0.5 1]);

% % add data from experiments
data=[0.757 0.673; %X2
    0.704 0.645]; %X3

for xi=1:2
    plot(data(xi,1), data(xi,2), 'ro');
end

set(gcf,'color','w','units','centimeters','pos',[5 5 10 10]);

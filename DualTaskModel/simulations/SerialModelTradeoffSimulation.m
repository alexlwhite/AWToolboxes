%% simulate accuracy tradeoffs using the all-or-none serial model
% Alex L. White

%Note, this script can take a long time to run, depending on nSimTrials

resFolder = fileparts(which('SerialModelTradeoffSimulation'));

nSimTrials = 1000000;

%% simulation of basic all-or-none serial model, for the prediction on graph

pProcessBoth = 0; %never process both stimuli
pStim1First = 0.5; %no bias to one side or the other

singleTaskDprimes = 0:0.2:4;
nDs = length(singleTaskDprimes);

accs = NaN(2, nDs);

for di = 1:nDs
    attnd1Mean = singleTaskDprimes(di);
    attnd2Mean = attnd1Mean;
    [~, ~, ~, ~, tradeoffAgs] = SerialDualTaskModel(attnd1Mean, attnd2Mean, pStim1First, pProcessBoth, nSimTrials);
    accs(:,di) = tradeoffAgs;
end

predT = table;
predT.d = singleTaskDprimes';
predT.AgGivenOtherSideIncorrect = accs(1,:)';
predT.AgGivenOtherSideCorrect = accs(2,:)';
predT.pProcessBoth = ones(nDs,1)*pProcessBoth;
predT.pStim1First = ones(nDs,1)*pStim1First;
predT.nSimTrials = ones(nDs,1)*nSimTrials;

save(fullfile(resFolder,'SimulatedAccuracyTradeoff_AllOrNoneSerial.mat'), 'predT');


%% plot
hues = [0.33 0.67];
sats = [0.8 0.8];
vals = [0.7 1];

colrs = hsv2rgb([hues' sats' vals']);

figure; hold on;
plot([0.5 1], [0.5 1],'k--');
plot(accs(1,:), accs(2,:), '.-', 'Color',colrs(1,:));


xlabel('Ag | other side incorrect');
ylabel('Ag | other side correct');

xlim([0.5 1]); ylim([0.5 1]);

for xi=1:2
    plot(data(xi,1), data(xi,2), 'ro');
end

set(gcf,'color','w','units','centimeters','pos',[5 5 10 10]);




%% Also, simulate the tradeoff pattern as a function of the bias to process one side or the other (pStim1First)

%pDualBoth: on what proportion of dual-task trials both sides can be processed
pDualBoth = 0;

singleTaskDprimes = 0:0.25:3.5;
nDs = length(singleTaskDprimes);

%Bias towards one side in dual-task (on trials when only 1 is attended)
%implemented by pStim1First:
%Ranges from 0 to 1, with 0 meaning never attend to the left, 0.5 being
%perfectly balanced, and 1 meaning always attend to the left
pStim1First = 0:0.25:1;
nBs = length(pStim1First);

accs = NaN(2, nDs, nBs);
for di = 1:nDs
    attnd1Mean = singleTaskDprimes(di);
    attnd2Mean = attnd1Mean;
    for bi = 1:nBs
        [~, ~, ~, ~, tradeoffAgs] = SerialDualTaskModel(attnd1Mean, attnd2Mean, pStim1First(bi), pDualBoth, nSimTrials);
        accs(:,di,bi) = tradeoffAgs;
    end
end

% %  plot
hues = linspace(0,1,nBs);
sats = ones(1,nBs)*0.8;
vals = linspace(0.3, 1, nBs);

colrs = hsv2rgb([hues' sats' vals']);

figure; hold on;
plot([0.5 1], [0.5 1],'k--');
hs = NaN(1,nBs);
legLabs = cell(1,nBs);
for bi=1:nBs
    hs(bi)=plot(accs(1,:,bi), accs(2,:,bi), '.-', 'Color',colrs(bi,:));
    legLabs{bi} = sprintf('%.2f', pStim1First(bi));
end

% % add data from experiments
data=[0.704 0.645; %X1
      0.757 0.673]; %X2
    
for xi=1:2
    plot(data(xi,1), data(xi,2), 'ko');
end

xlim([0.5 1]); ylim([0.5 1]);

xlabel('Ag | other side incorrect');
ylabel('Ag | other side correct');

legend(hs,legLabs,'Location','NorthWest');

title('Varying bias to process side 1');
set(gcf,'color','w','units','centimeters','pos',[5 5 10 10]);

% Conculsion: the tradeoff pattern is not affected by the bias to process
% one side or the other

%% Now simulate the same thing, but varying the probability of processing both stimuli on dual-task trials (pDualBoth)
singleTaskDprimes = 0:0.25:3.5;
nDs = length(singleTaskDprimes);

%probability of processing both words on dual-task trials:
pBoths = 0:0.2:1;
nBs = length(pBoths);

pStim1First = 0.5; %no bias

accs = NaN(2, nDs, nBs);

for di = 1:nDs
    attnd1Mean  = singleTaskDprimes(di);
    attnd2Mean  = attnd1Mean;
    for bi = 1:nBs
        [~, ~, ~, ~, tradeoffAgs] = SerialDualTaskModel(attnd1Mean, attnd2Mean, pStim1First, pBoths(bi),nSimTrials);
        accs(:,di,bi) = tradeoffAgs;
    end
end

% % plot
hues = linspace(0,1,nBs);
sats = ones(1,nBs)*0.8;
vals = linspace(0.3, 1, nBs);

colrs = hsv2rgb([hues' sats' vals']);

figure; hold on;
plot([0.5 1], [0.5 1],'k--');
hs = NaN(1,nBs);
legLabs = cell(1,nBs);
for bi=1:nBs
    hs(bi)=plot(accs(1,:,bi), accs(2,:,bi), '.-', 'Color',colrs(bi,:));
    legLabs{bi} = sprintf('%.2f', pBoths(bi));
end

xlabel('Ag | other side incorrect');
ylabel('Ag | other side correct');

xlim([0.5 1]); ylim([0.5 1]);

legend(hs,legLabs,'Location','NorthWest');
title('Varying probability of processing both');

for xi=1:2
    plot(data(xi,1), data(xi,2), 'ro');
end

set(gcf,'color','w','units','centimeters','pos',[5 5 10 10]);

% Conculsion: the tradeoff effect is diminished the more often both
% stimuli get processed
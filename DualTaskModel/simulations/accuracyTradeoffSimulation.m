%% simulate accuracy tradeoffs, aka contingent Ag

resFolder = '/Users/alexwhite/Dropbox/PROJECTS/DualTaskModel/simulationResults';

nSimTrials = 50000;
modelType = 1; %1=serial, 2=parallel

%Extra capacity parameter: on what proportion of dual-task trials both
%sides can be processed
pDualBoth = 0;

%Congruency effects implemented by pFlip2Side2:
%On some proportion of trials, attend to side 2 only, but without
%realizing it, so report that side even if asked about the right
%pFlip2Side2 ranges from -.5 (flip to side 1 half the time) to 0.5 (flip to side 2)
pFlip2Side2 = 0;

singleTaskDprimes = 0:0.25:3.5;
nDs = length(singleTaskDprimes);

%Bias towards one side in dual-task (on trials when only 1 is attended)
%implemented by dualTaskAttnBias:
%Ranges from 0 to 1, with 0 meaning never attend to the left, 0.5 being
%perfectly balanced, and 1 meaning always attend to the left
dualTaskAttnBiases = 0:0.25:1;
nBs = length(dualTaskAttnBiases);

accs = NaN(2, nDs, nBs);

for di = 1:nDs
    attnd1Mean = singleTaskDprimes(di);
    attnd2Mean = attnd1Mean;
    for bi = 1:nBs        
        [~, ~, ~, ~, contingentAgs, ~] = GeneralDualTaskModel(modelType,attnd1Mean, attnd2Mean, pFlip2Side2, dualTaskAttnBiases(bi), pDualBoth, nSimTrials);
        accs(:,di,bi) = contingentAgs(1,1,:);
    end
end

% interesting note: the bias to report one side or the other doesn't change
% the size of the tradeoff 

save(fullfile(resFolder,'SimulatedContingentAcc_dprimeAndBias.mat'), 'accs','singleTaskDprimes','dualTaskAttnBiases');
%% plot
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
  legLabs{bi} = sprintf('%.2f', dualTaskAttnBiases(bi));
end

% % add data from experiments
data=[0.757 0.673; %new X2 
     0.704 0.645]; %new X3
 
 for xi=1:2
     plot(data(xi,1), data(xi,2), 'ko');
 end
 
xlim([0.5 1]); ylim([0.5 1]); 

 
xlabel('Ag | other side incorrect');
ylabel('Ag | other side correct');

legend(hs,legLabs,'Location','NorthWest');

title('Varying bias to process side 1');



set(gcf,'color','w','units','centimeters','pos',[5 5 10 10]);
figTitle = fullfile(resFolder, 'AccTradeoff_VaryingAttnSideBias.eps');
exportfig(gcf,figTitle,'Format','eps','bounds','loose','color','rgb','LockAxes',0,'FontMode','fixed','FontSize',12);

 
 %% do it all again, but vary pDualBoth
 nSimTrials = 50000;

modelType = 1; %1=serial, 2=parallel
pFlip2Side2 = 0;

singleTaskDprimes = 0:0.25:3.5;
nDs = length(singleTaskDprimes);

%probability of processing both words on dual-task trials:
 pBoths = 0:0.2:1; 
 nBs = length(pBoths);
 
 dualTaskAttnBias = 0.5; %no bias
 
 accs = NaN(2, nDs, nBs);

for di = 1:nDs
    attnd1Mean  = singleTaskDprimes(di);
    attnd2Mean  = attnd1Mean;
    for bi = 1:nBs        
        [~, ~, ~, ~, contingentAgs, ~] = GeneralDualTaskModel(modelType,attnd1Mean, attnd2Mean, pFlip2Side2, dualTaskAttnBias, pBoths(bi),nSimTrials);
        accs(:,di,bi) = contingentAgs(1,1,:);
    end
end

save(fullfile(resFolder,'SimulatedContingentAcc_dprimeAndPBoth.mat'), 'accs','singleTaskDprimes','dualTaskAttnBiases');
%% plot
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

% % add data from experiments
data=[0.757 0.673; %new X2 
     0.704 0.645]; %new X3
 
 for xi=1:2
     plot(data(xi,1), data(xi,2), 'ro');
 end
 
set(gcf,'color','w','units','centimeters','pos',[5 5 10 10]);
figTitle = fullfile(resFolder, 'AccTradeoff_VaryingProbProcessBoth.eps');
exportfig(gcf,figTitle,'Format','eps','bounds','loose','color','rgb','LockAxes',0,'FontMode','fixed','FontSize',12);

%% basic serial and fixed capacity model comparison 

nSimTrials = 1000000;
modelTypes = 1:2; %1=serial, 2=parallel
modelLabels = {'All-or-none serial','Fixed cap. parallel'};

nModels = length(modelTypes);

%Extra capacity parameter: on what proportion of dual-task trials both
%sides can be processed (for serial model); 
%or how much extra capacity is there in dual-task trials (for parallel
%model)
extraCapacity = 0; 

%Congruency effects implemented by pFlip2Side2:
%On some proportion of trials, attend to side 2 only, but without
%realizing it, so report that side even if asked about the right
%pFlip2Side2 ranges from -.5 (flip to side 1 half the time) to 0.5 (flip to side 2)
pFlip2Side2 = 0;

singleTaskDprimes = 0:0.25:3.5;
nDs = length(singleTaskDprimes);




%Bias towards one side in dual-task (on trials when only 1 is attended)
%implemented by dualTaskAttnBias:
%Ranges from 0 to 1, with 0 meaning never attend to the left, 0.5 being
%perfectly balanced, and 1 meaning always attend to the left
dualTaskAttnBias = 0.5;

accs = NaN(2, nDs, nModels);

for di = 1:nDs
    attnd1Mean = singleTaskDprimes(di);
    attnd2Mean = attnd1Mean;
    for mi = 1:nModels       
        [~, ~, ~, ~, contingentAgs, ~] = GeneralDualTaskModel(modelTypes(mi),attnd1Mean, attnd2Mean, pFlip2Side2, dualTaskAttnBias, extraCapacity, nSimTrials);
        accs(:,di,mi) = contingentAgs(1,1,:);
    end
end


save(fullfile(resFolder,'SimulatedContingentAcc_SerialVsFixedCapParallel.mat'), 'accs','singleTaskDprimes','modelTypes');



% % plot
hues = [0.33 0.67]; 
sats = [0.8 0.8]; 
vals = [0.7 1];

colrs = hsv2rgb([hues' sats' vals']);

figure; hold on; 
plot([0.5 1], [0.5 1],'k--');
hs = NaN(1,nModels);
for mi=1:nModels
   hs(mi)=plot(accs(1,:,mi), accs(2,:,mi), '.-', 'Color',colrs(mi,:));
end

xlabel('Ag | other side incorrect');
ylabel('Ag | other side correct');

xlim([0.5 1]); ylim([0.5 1]); 

legend(hs,modelLabels,'Location','NorthWest');
 title('Serial vs parallel');

% % add data from experiments
data=[0.757 0.673; %new X2 
     0.704 0.645]; %new X3
 
 for xi=1:2
     plot(data(xi,1), data(xi,2), 'ro');
 end
 
set(gcf,'color','w','units','centimeters','pos',[5 5 10 10]);
figTitle = fullfile(resFolder, 'AccTradeoff_Serial vs parallel.eps');
exportfig(gcf,figTitle,'Format','eps','bounds','loose','color','rgb','LockAxes',0,'FontMode','fixed','FontSize',12);

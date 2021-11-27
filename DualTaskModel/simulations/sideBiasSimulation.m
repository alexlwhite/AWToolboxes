%% simulate division of resources between the two sides 

nSimTrials = 50000;
modelType = 2; %1=serial, 2=parallel

%sides can be processed
pDualBoth = 0;

%Congruency effects implemented by pFlip2Side2:
%On some proportion of trials, attend to side 2 only, but without
%realizing it, so report that side even if asked about the right
%pFlip2Side2 ranges from -.5 (flip to side 1 half the time) to 0.5 (flip to side 2)
pFlip2Side2 = 0;


singleTask1Dprime = AgToDprime(0.8);
singleTask2Dprime = AgToDprime(0.85);

%Bias towards one side in dual-task (on trials when only 1 is attended)
%implemented by dualTaskAttnBias:
%Ranges from 0 to 1, with 0 meaning never attend to the left, 0.5 being
%perfectly balanced, and 1 meaning always attend to the left
dualTaskAttnBiases = 0:0.1:1;
nBs = length(dualTaskAttnBiases);

accs = NaN(2,2,nBs);

congI = 1; %ignore congruency for now
for bi = 1:nBs
    [Ags, valsByIndex] = GeneralDualTaskModel(modelType,singleTask1Dprime, singleTask2Dprime, pFlip2Side2, dualTaskAttnBiases(bi), pDualBoth, nSimTrials);
    %extract data
    for condSide=1:3 %single-task 1, single-task 2, dual-task
        if condSide<3 %single-task
            sideI = condSide;
            thisCond = 1;
        else %dual-task
            sideI=1:2;
            thisCond = 2;
        end
        accs(thisCond, sideI, bi) = squeeze(Ags(condSide, sideI, congI));
    end
end

%% plot 

hues = linspace(0,1,nBs); 
sats = ones(1,nBs)*0.8; 
vals = linspace(0.3, 1, nBs);

colrs = hsv2rgb([hues' sats' vals']);


figure; hold on; 

%plot predicted AOC to compare 
as = mean(accs,3);
as(2,:) = accs(2,:,ceil(nBs/2));

es = NaN(2,2);
plotOpt.edgeColors = zeros(2,3); 
plotOpt.fillColors = plotOpt.edgeColors; 
plotOpt.fillColors(2,:) = 1; 
plotOpt.doLegend = false;
plotOpt.sideLabels = {'stim1','stim2'};
plotOpt.plotSerialPrediction = false;

plotAOCWithPredictions(as, es, plotOpt);

for bi=1:nBs
   plot(0.5, accs(1,1, bi), 'o','MarkerEdgeColor',colrs(bi,:),'MarkerFaceColor', colrs(bi,:)); 
   plot(accs(1,2,bi), 0.5, 'o','MarkerEdgeColor',colrs(bi,:),'MarkerFaceColor', colrs(bi,:));  
   plot(accs(2,2,bi), accs(2,1,bi), 'o','MarkerEdgeColor',colrs(bi,:),'MarkerFaceColor', 'w'); 
    
end

xlim([0.5 1]); ylim([0.5 1]); 

title('Varying the bias to attend to one side');

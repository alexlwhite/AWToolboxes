%% This script creates an imaginary psychometric function in a yes/no visual 
%% detection ask, and simulates a SIAM staircase procedure to estimate the treshold 
% 
% by Alex L. White, 2012-2018, at the University of Washington 

clear; 
home;


nSim = 1000; %number of simulations 

nTrls = 200; %total number of trials per simulation 

%% Determine which contrast levels are available

scr.black = 0;
scr.white = 255;
scr.bgColor = scr.black+round(0.5*(scr.white-scr.black));  % background color

%Which contrasts are really available
cUps=(scr.bgColor+1):scr.white;
cDns=(scr.bgColor-1):-1:scr.black; 

if length(cDns)>length(cUps)
    cDns=cDns(1:length(cUps));
end

contrastSet=(cUps-cDns)./(cUps+cDns);

%or just set it:
%contrastSet = logspace(log10(0.01), log10(1),100);

%% Set up the 'true' psychometric function as a Naka-Rushton: d' as a function of contrast

Rmax = 3; %upper asymptote
C50  = 0.1; %C50 is the threhshold: the contrast level where d' reaches half the upper asymptote
n    = 2; %slope

nrParams = [Rmax C50 n];

cx  = logspace(log10(0.001),log10(1),6);
cxl = log10(cx);

psi = nakaRushton(nrParams,contrastSet);

figure(1); clf;
subplot(2,1,1),
set(gcf,'Position',[5 5 500 850]);

plot(log10(contrastSet),psi);
hold on; 
%plot c50
c50h = plot(log10([C50 C50]), [0 nakaRushton(nrParams, C50)],'g-','LineWidth',3);


xlim([log10(contrastSet(1)/2) log10(contrastSet(end)*1.5)]); 
ylim([0 Rmax*1.1]);
set(gca,'XTick',cxl);

for xi=1:length(cx)
    contrLabs{xi} = sprintf('%.3f',cx(xi));
end
set(gca,'XTickLabel',contrLabs);
xlabel('Contrast');
ylabel('d''');

title('Psychometric function and staircase estimate','FontSize',15);
set(gca,'FontSize',12);

%% Set up SIAM staricase 
%total number of reversals before stopping staircase
goalReversals = 20; 

%starting step size
startStep = log10(contrastSet(2))-log10(contrastSet(1));
%the final step size should be bigger than the biggest difference between
%availabe contrasts, to avoid getting stuck!  

%other paramters 
goalAccuracy  = 0.5; %maximum reduced hit rate, ranges from 0-1
threshType    = 1; %mean
revsToHalveI  = [1 2];
revsToReset   = 20;
nStuckToReset = 5;


%% Simulate using staircase 
threshs  = zeros(1,nSim); 
totalNTs = zeros(1,nSim); 

for simi = 1:nSim
    %set the starting level, with a noisy guess of true threshold: 
    %the true c50 + noise
    startC = C50+rand*.4;
    
    %initialize the staircase
    ss = initSIAM(goalAccuracy,startStep,log10(startC),log10(contrastSet),revsToHalveI,revsToReset,nStuckToReset);
    %for ti=1:nTrls
    %keep going through trials until reach the required number of reversals
    while ss.nreversals<goalReversals
        %whether a target is present 
        pres=CoinFlip(1,0.5); 
        
        %which contrast level the staircase recommends. 
        %note that the staircase recommends a log level, so here we un-log it
        contrast=10^ss.intensity; 
        
        %the true d' level for this contrast level
        d=nakaRushton(nrParams,contrast); 
        
        %how to get from a d' level to hit and false alarm rates, assuming neutral
        %criterion?
        %Neutral criterion means that hit rate = 1 - false alarm rate
        %so d=norminv(HR)-norminv(FAR)
        %   d=2*norminv(HR); 
        %   d/2=norminv(HR); 
        %   HR = normcdf(d/2); 
        
        %determine the observer's response (target present vs absent)
        hitRate=normcdf(d/2); 
        if pres
            resp=CoinFlip(1,hitRate); 
        else
            resp=CoinFlip(1,1-hitRate); 
        end
        
        ss = updateSIAM(ss,pres,resp); 
    end
    
    %estimate threshold as median of all intensities after the first five reversals
    threshs(simi) = estimateSIAM(ss,threshType);

    totalNTs(simi)= ss.tnum; 
    
    %plot the first staircase:
    if simi==1
        subplot(2,1,2); hold on;
        plotSIAM(ss,threshType); 
        title('Example staircase');
    end
end

%% Plot the mean threshold estimate on the psychometric function 
meanThresh = mean(threshs); 

%and the 95% confidence inteval 
threshCICs=[quantile(threshs,0.025) meanThresh quantile(threshs,0.975)];

figure(1); subplot(2,1,1); hold on;
threshDs = nakaRushton(nrParams,10.^threshCICs);
lineWidths=[1 2 1];
for cdi=1:3
    threshH(cdi)=plot([threshCICs(cdi) threshCICs(cdi)],[0 threshDs(cdi)],'k-','LineWidth',lineWidths(cdi));
end

legend([c50h threshH(2)],'c50','Mean Threshold + 95% CI','Location','NorthWest');




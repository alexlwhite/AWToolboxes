%first, make a naka-rushton function describing d' as a function of
%contrast 
%close all; 
clear; 
home;


nsim=1;
nts=200;

%% Contrasts

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

Rmax = 3;
C50 = 0.1; 
n = 2;

nrParams=[Rmax C50 n];

%contrastSet = logspace(log10(0.01), log10(1),100);

cx = logspace(log10(0.001),log10(1),6);
cxl = log10(cx);

psi = nakaRushton(nrParams,contrastSet);

figure; 
if nsim==1
    subplot(2,1,1), 
    set(gcf,'Position',[5 5 500 850]);
else
    set(gcf,'Position',[5 5 600 600]);
end
plot(log10(contrastSet),psi);
hold on; 
%plot c50
c50h=plot(log10([C50 C50]), [0 nakaRushton(nrParams, C50)],'g-','LineWidth',3);


xlim([log10(contrastSet(1)/2) log10(contrastSet(end)*1.5)]); 
ylim([0 Rmax*1.1]);
set(gca,'XTick',cxl);

for xi=1:length(cx)
    contrLabs{xi}=sprintf('%.3f',cx(xi));
end
set(gca,'XTickLabel',contrLabs);
title('Psychometric function and staircase estimate','FontSize',15);
set(gca,'FontSize',12);



goalReversals=20; %total number of reversals to keep doing

startStep=log10(contrastSet(2))-log10(contrastSet(1));
%startStep=log10(0.1)-log10(0.05);
%the final step size should be bigger than the biggest difference between
%availabe contrasts, to avoid getting stuck!  


threshs=zeros(1,nsim); 
totalNTs=zeros(1,nsim); 

revsToHalveI = [1 2];
revsToReset = 20;

for simi = 1:nsim
    startC = C50*(0.8+rand*.4); %starting level is c50 with 20% error
    
    %initialize staircase
    ss = initSIAM(0.5,startStep,log10(startC),log10(contrastSet),revsToHalveI,revsToReset);
    for ti=1:nts
    %while ss.nreversals<goalReversals
        pres=CoinFlip(1,0.5); 
        contrast=10^ss.intensity; 
        d=nakaRushton(nrParams,contrast); 
        
        %how to get from a d' level to hit and false alarm rates, assuming neutral
        %criterion?
        %neutral criterion means that hit rate = 1 - false alarm rate
        %so d=norminv(HR)-norminv(FAR)
        %   d=2*norminv(HR); 
        %   d/2=norminv(HR); 
        %   HR = normcdf(d/2); 
        
        hitRate=normcdf(d/2); 
        if pres
            resp=CoinFlip(1,hitRate); 
        else
            resp=CoinFlip(1,1-hitRate); 
        end
        
        ss = updateSIAM(ss,pres,resp); 
    end
    
    %estimate threshold as median of all intensities after the first five
    %reversals
    threshs(simi) = estimateSIAM(ss,1);

    totalNTs(simi)= ss.tnum; 
end

%show the mean threshold estimate on the psychometric function 
meanThresh=mean(threshs); 
stdThresh=std(threshs); 

threshCICs=[quantile(threshs,0.025) meanThresh quantile(threshs,0.975)];


%Plot 95% confidence interval
figure(1); hold on; 
threshDs = nakaRushton(nrParams,10.^threshCICs);
lineWidths=[1 2 1];
for cdi=1:3
    threshH(cdi)=plot([threshCICs(cdi) threshCICs(cdi)],[0 threshDs(cdi)],'k-','LineWidth',lineWidths(cdi));
end

legend([c50h threshH(2)],'c50','Mean Threshold','Location','NorthWest');

xlabel('Contrast');
ylabel('d''');

if nsim==1
    %plot staircase 
    
    figure(1); subplot(2,1,2); hold on;
        
    is=10.^ss.ints;

    plot(1:ss.tnum,is,'b-'); 
    
    %plot hits
    hitTrls = find(ss.pres & ss.resp);
    htH=plot(hitTrls, is(hitTrls), 'g.','MarkerSize',10);
    
    missTrls = find(ss.pres & ~ss.resp);
    msH=plot(missTrls, is(missTrls), 'r.','MarkerSize',10);
    
    crTrls = find(~ss.pres & ~ss.resp);
    crH=plot(crTrls, is(crTrls), 'b.','MarkerSize',10);
    
    faTrls = find(~ss.pres & ss.resp);
    faH=plot(faTrls, is(faTrls), 'y.','MarkerSize',10);
    
    for ri=ss.reversalTs
        plot([ri ri],[0 is(ri)],'k-');
    end
    for rsti=ss.resetT
        plot([rsti rsti],[0 is(rsti)],'-','Color',[.8 0 .7],'LineWidth',1.5);
    end
    
    plot([0 ss.tnum],10.^[meanThresh meanThresh],'r-');
    
    title('Staircase','FontSize',15);
    xlabel('Trial'); 
    ylabel('Contrast'); 
    
    legend([htH msH crH faH],'Hit','Miss','Corr Reject','False Alarm','Location','NorthEast');
    
    set(gca,'FontSize',12);
    
end

mean(totalNTs)

%test the function withinSubjectSEMs on a fake data set 

means = [5 6 9 7]; %initial condition means controlling effect sizes
ncond = length(means);
nsubj = 10;
%generate data with no systematic across-subject differences: 
dat1 = randn(nsubj,ncond)+repmat(means,10,1); 

%generate overall across-subject differences, as offsets applied to all of each subject's data points
subjMeans = 5*randn(nsubj,1); 

%generate full data with consistent within-subject differences across conditions but big across-subject mean differences as well
dat2 = dat1+repmat(subjMeans,1,ncond);  

meanDat = mean(dat2,1); %mean data 
sems = std(dat2,0,1)/sqrt(nsubj); %standard errors 

%compute within-subject standard errors following Morey (2008) 
withinSEMs = withinSubjectSEMs(dat2); 

%% plot mean data with 2 types of error bars
figure; hold on;
for ci=1:ncond
    h1=plot([ci ci],meanDat(ci)+[-1 1]*sems(ci),'b-','LineWidth',3); 
    h2=plot([ci ci],meanDat(ci)+[-1 1]*withinSEMs(ci),'r-','LineWidth',2);
    h3=plot(ci,meanDat(ci),'k.','MarkerSize',20);
end

xlim([0 ncond+1]); ylim([0 ceil(0.5+max(meanDat+sems))])
xlabel('Condition'); 
ylabel('Measurement'); 
legend([h3 h1 h2],{'mean','+/- 1 SEM','+/- 1 within-subject SEM'},'Location','NorthWest');
legend boxoff;
axis square;


% function [A, HRs, FARs] = ROC(p, s)
% Alex White, April 2017
% 
% ROC analysis: determine area under the ROC curve in a signal-detection
% scenario. 
% 
%Inputs
% - p: 1xn vector of strengths of evidence for signal
%      present events on each of n trials
% - s: 1xn vector of Boolean variables, for each "trial" whether there was
%      a signal present
% 
%Outputs 
% - A: Area under the ROC curve 
% - HRs: 1xC vector of hit rates for each of C criteria, which linearly
% decrease from just above the highest p to below just the lowest value of
% p. That means HRs monotonically increase. 
% - FARs: 1xC vector of false alarm rates. 

function [A, HRs, FARs] = ROC(p, s)

if length(unique(s))<=1
    A = NaN;
    fprintf(1,'\n(ROC): Area under ROC is undefined because either all S or all ~S\n');
elseif length(unique(round(p*10000000)/10000000))<=1 
    A = NaN;
    fprintf(1,'\n(ROC): Area under ROC undefined because there is only one value of p\n');
else
    
%Determine criteria c
pmin = min(p); 
pmax = max(p);


%number of criteria
nc = length(p); 
nc(nc<4)=4;

%distance between criteria: 
cw = (pmax-pmin)/(nc-3);

%Construct vectors cs of criteria, in decreasing  order so HR and FAR 
%monotonically increase 
%Also, Add one c above the max p, and one below the min p, so that HR and FAR
%both to to 0 and 1 at the extremes 
cs=linspace((pmax+cw),(pmin-cw), nc);

%Compute hit and false alarm rates for each c
HRs = zeros(1,nc);
FARs = zeros(1,nc); 

for ci = 1:nc
    HRs(ci) = mean(p(s==1)>cs(ci));
    FARs(ci) = mean(p(s==0)>cs(ci)); 
end


%plotting: 
% figure; subplot(1,2,1); hold on;
% plot(cs,HRs,'r.-');
% plot(cs,FARs,'g.-');
% xlabel('Criterion');
% ylabel('HR or FAR');
% axis square;
% 
% subplot(1,2,2);
% hold on;
% plot([0 1],[0 1],'k-')
% plot(FARs,HRs,'.-');
% xlabel('FAR'); ylabel('HR');
% axis square;

%Compute area under ROC curve: 
A = computeAROC(HRs,FARs);
end
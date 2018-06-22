%% function [contTable,k,corrRho,corrP] = computeEventContingency(xs)
% Computes contingency tables and correlation between two sets of binary
% event types 
% Inputs: 
% xs: a Nx2 matrix of event outcomes, 0 or 1. Each row is one time-point or
% 'trial'. Column 1 = success or failure of event type 1; Column 2 = success or failure of event type 2. 
% 
% Outputs: 
% - contTable: contingency table for the two event types. Format: 
%     [p(success1 & success2) p(success1 & fail2); 
%      p(fail1 & success2)    p(fail1 & fail2)]
% - k: the amount by which p(success1 & success2) differs from that
% predicted by independence. That should be the same for p(fail1 & fail2),
% and -1* the difference for the other two joint probabilities. 
% -corrRho: correlation coefficient for correlation of the 2 columns in xs
% -corrP: the pval for corrRho

function [contTable,k,corrRho,corrP] = computeEventContingency(xs)

%marginal probabilities of success for event 1 (left column of xs) and
%event 2 (right column of xs)
marginalPs = mean(xs,1); 

%fill in contingency table if event 1 and event 2 were independent 
indepTable = NaN(2,2); 

indepTable(1,1) = marginalPs(1)*marginalPs(2);  %both success 
indepTable(1,2) = marginalPs(1)*(1-marginalPs(2)); %left success, right fail
indepTable(2,1) = (1-marginalPs(1))*(marginalPs(2)); %left fail, right success
indepTable(2,2) = (1-marginalPs(1))*(1-marginalPs(2)); %both fail 

%compute actual contingency table 
contTable = NaN(2,2); 
contTable(1,1) = mean(xs(:,1)==1 & xs(:,2)==1); %both sides success 
contTable(1,2) = mean(xs(:,1)==1 & xs(:,2)==0); %left success, right fail 
contTable(2,1) = mean(xs(:,1)==0 & xs(:,2)==1); %left fail, right success 
contTable(2,2) = mean(xs(:,1)==0 & xs(:,2)==0); %both sides fail 

%compute deviation from independence: 
ks = contTable - indepTable;
k = ks(1,1);

%also compute linear correlation 
[corrRho, corrP] = corr(xs(:,1),xs(:,2));
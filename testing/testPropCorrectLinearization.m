% % generate d' values 
ds = 0:0.1:5;

% % convert to proportion correct (derived from basic signal detection
% theory, assuming a neutral criterion, such that hit rate = correct reject rate)
ps = normcdf(ds/sqrt(2));

% % plot d' vs p(correct)
figure; subplot(1,3,1); hold on;
plot(ds, ps, 'b-'); 
xlabel('d''');
ylabel('p(correct)'); 
xlim([0 5]); ylim([0.5 1]);
title('d'' vs. prop correct');

% % so, d' is non-linearly related to p(correct). In psychophysics, d' is the
% what we want to really measure: the SNR of the information the brain uses
% to make a decision. Here we see that p(correct) gets compressed at the
% upper range. 

% % So can we fix that by z-scoring the proportions? 
% first try with matlab's z-score function. That basically returns each
% number in a vector transformed to reflect the number of SDs away from the
% mean of the vector. 
psZ = zscore(ps);

subplot(1,3,2); hold on;
plot(ds, psZ, 'r-');
xlabel('d''');
ylabel('zscore(pCorrect)'); 
title('d'' vs z-scored pcorrect');

% % but sometimes z-transforming just means using the inverse of the normal
% CDF (norminv). That would partly invert the equation used to go from d'
% to p(correct) and should linearize. 
psZ2 = norminv(ps); 

subplot(1,3,3); hold on;
plot(ds, psZ2, 'g-');
xlabel('d''');
ylabel('norminv(pCorrect)'); 
title('d'' vs norminv(pcorrect)');

%yes, transforming p(correct) with norminv returns value that are linear
%with respect to d'. 

%% and what about proportion errors, rather than proportion correct? 
pE = 1-ps; %proportion error 
pEZ = norminv(pE); %z-transformed

figure;
plot(fliplr(ds), pEZ);

xlabel('d'''); 
ylabel('norminv(prop error)');

%this seems to work for proportion errors too... they're linear with
%respect to corresponding values of d'. 
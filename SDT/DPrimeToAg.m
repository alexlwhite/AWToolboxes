%% function Ag = DPrimeToAg(dprime)
% Converts between dprime and area under the ROC curve (Ag). 
% Treats Ag like unbiased p(correct) in a 2AFC task. 
function Ag = DPrimeToAg(dprime)

Ag = normcdf(dprime/sqrt(2));

%% What's the relationship betwen dprime and area under the ROC curve? 
%First let's test out computing area under the curve, Ag, for a range of d'
% %values. 
% % 
% %classic SDT model, assuming a value of d', the separation of two distributions with unit variance 
% %values of dprime: 
% ds =  0:0.1:7; 
% 
% %then we vary criterion and compute hit and false alarm rates along the way
% cs = 10:(-0.01):-10; 
% 
% %false alarm rates (don't depend on d')
% frs = 1-normcdf(cs,0,1); 
% 
% %Compute Ag corresponding to each value of d'
% Ags = NaN(1,length(ds));
% PCs = Ags;
% for di = 1:length(ds)
%     %hit rates given this value of dprime   
%     hrs = 1-normcdf(cs,ds(di),1);
%     %the compute area under the curve by summing up area of the trapezoids:
%     Ags(di) = computeAROC(hrs,frs);
%     %proportion correct with neutral criterion
%     neutralHR = 1-normcdf(ds(di)/2, ds(di), 1); %hit ratewhen c is midway between 0 and d'
%     neutralCR = normcdf(ds(di)/2, 0, 1); %correct reject rate 
%     PCs(di) = mean([neutralHR neutralCR]); %proportion correct is mean of HR and CR
% end
% 
% %Now, from basic SDT math, we have a way to convert from p(corr) to
% %d-prime, assuming neutral criterion. 
% %The formula is:
% %   d' = sqrt(2)*z(pc)
% % and  applies for 2-interval forced-choice experiments. The sqrt(2) factor basically corrects for the fact that 
% % the observer has 2x as much information. z(x) is the same as matlab's
% % norminv(x), inverse of normCDF. (But why do we use the formula for 2IFC here? Perhaps because 1 stim can give evidence for either option)
% % 
% %Can we treat Ag (area under curve) just like 2IFC p(corr), to translate from 
% % from Ag to dprime? 
% predDs = sqrt(2)*norminv(Ags); %convert to d' (assuming Ag is like p(corr))
% predAgs = normcdf(ds/sqrt(2)); %convert dprime to Ag also
% 
% predPC = normcdf(ds/2);
% 
% figure;  hold on;
% plot(ds,Ags,'b-','LineWidth',4);
% plot(predDs,predAgs,'g-','LineWidth',2);
% plot(ds, PCs, 'r-', 'LineWidth', 4);
% plot(ds, predPC, 'y-', 'LineWidth', 1.5);
% legend({'true Ag', 'predicted Ag','true PC','predicted PC'},'Location', 'NorthWest');
% xlabel('d-prime'); 
% ylabel('Ag');

%So, yes, we can summarize the relationship between Ag (area under ROC
%curve) and dprime as follows: 
%% dprime = sqrt(2)*norminv(Ag)
%% Ag = normcdf(dprime/sqrt(2))

%% But note that the relationship between dprime and a one-interval,
%2-alternative forced-choice proportion correct is different! (assuming
%unbiased criterion)
%dprime = 2*norminv(PC);
%PC = normcdf(dprime/2);
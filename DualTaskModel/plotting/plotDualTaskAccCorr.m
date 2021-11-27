function plotDualTaskAccCorr(rhos,condLabs)

figure; hold on;
xs = 1:length(rhos);
plot([0 max(xs)],[0 0],'k-');
plot(xs, rhos, 'b.','MarkerSize',25);

set(gca,'XTick',xs,'XTickLabel',condLabs); 
ylabel('Corr. coeff.'); 

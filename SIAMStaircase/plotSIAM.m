function [] = plotSIAM(ss,threshType)

hold on;

is=10.^ss.ints;
thresh = estimateSIAM(ss,threshType);

plot(1:ss.tnum,is,'b-');

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

for sri = ss.revStableStepsTs
    plot([sri sri],[0 is(sri)],'k-','LineWidth',2);
end

for rsti=ss.resetT
    plot([rsti rsti],[0 is(rsti)],'-','Color',[.8 0 .7],'LineWidth',2);
end

plot([0 ss.tnum],10.^[thresh thresh],'r-');

%title('Staircase','FontSize',15);
xlabel('Trial');
ylabel('Intensity');

legend([htH msH crH faH],'Hit','Miss','Corr Reject','False Alarm','Location','NorthEast');

set(gca,'FontSize',12);
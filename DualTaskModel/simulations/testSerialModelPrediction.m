%% test the simple serial model as described in the manuscript 
% Plot the diagonal line prediction in the AOC 

%proportion of trials both stimuli are processed
w = 0; 

%example accuracy levels for top and bottom stimulus 
singleTaskAgs = [0.82 0.87]; 

%vary the proportion of trials that stimulus 1 is processed 'first' 
vs = 0:0.1:1;
nPs = length(vs);

dualAgs = NaN(2,nPs);

for pi = 1:nPs
    
    for si=1:2 %task side
        if si==1
            pProcessFirst = vs(pi);
        else
            pProcessFirst = 1-vs(pi);
        end
        dualAgs(si,pi) = w*singleTaskAgs(si) + (1-w)*(pProcessFirst*singleTaskAgs(si) + 0.5*(1-pProcessFirst));
    end
end


figure; hold on;
plot(0.5, singleTaskAgs(1), 'o','MarkerSize',12,'MarkerEdgeColor','b','MarkerFaceColor','b');
plot(singleTaskAgs(2), 0.5, 'o','MarkerSize',12,'MarkerEdgeColor','b','MarkerFaceColor','b');

plot(dualAgs(2,:),dualAgs(1,:), 'k-');

xlim([0.5 1]); ylim([0.5 1]);

xlabel('Bottom accuracy'); 
ylabel('Top accuracy');
axis square

%% equivalently
% for pi = 1:nPs
%     
%     for si=1:2 %task side
%         if si==1
%             pProcess = (1-w)*vs(pi) + w;
%         else
%             pProcess = w+(1-w)*(1-vs(pi));
%         end
%         dualAgs(si,pi) = singleTaskAgs(si)*pProcess + 0.5*(1-pProcess);
%     end
% end
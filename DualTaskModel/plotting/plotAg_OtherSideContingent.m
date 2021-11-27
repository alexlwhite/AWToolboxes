function plotAg_OtherSideContingent(Ag,valsByIndex)

%colors 
hues = [0.4 0.4]; 
sats = [0.9 0.5];
vals = [0.5 1]; 

otherCorrCols = hsv2rgb([hues' sats' vals']);

datLineWidth = 3;
axLineWidth = 0.5;

%Y axis
ylims = [0.5 1];
yticks = [0.5:0.1:1];
ytickLabs = cell(1,length(yticks));
for yti=1:length(yticks)
    v = yticks(yti);
    if mod(yti,2)==1
        if (round(v)-v)==0
            frmt = '%.1f';
        elseif (round(10*v)-10*v) == 0
            frmt = '%.1f';
        elseif (round(100*v)-100*v) == 0
            frmt = '%.2f';
        else
            frmt = '%.3f';
        end
        ytickLabs{yti} = sprintf(frmt,v);
    else
        ytickLabs{yti} = '';
    end
end


bwid=.7;
bsep=.3;

figure; hold on;
task = 1;

sidesToPlot = 1:2;
otherSidePresToPlot = [0 1];
otherSideCorrToPlot = [0 1];

sideLabels = cell(1,2);
for si = 1:length(sidesToPlot)
    
    sideI = find(valsByIndex.sides==sidesToPlot(si));
    
    if sidesToPlot(si) == 1
        sideLabels{si} = 'Left';
    elseif sidesToPlot(si) == 2
        sideLabels{si} = 'Right';
    end
    
    ds = squeeze(Ag(sideI,:,:));
        
    subplot(1,2,si);
    hold on;
    
    bi = 0;
    bctr = -bsep;
    ctrs = [];
    groupctrs = [];
    
    for pi = 1:length(otherSidePresToPlot)
        otherPresI = find(valsByIndex.otherSidesPres==otherSidePresToPlot(pi));
        
        if otherSidePresToPlot(pi)==0
            otherPresLabs{pi} = 'Absent';
        elseif otherSidePresToPlot(pi)==1
            otherPresLabs{pi} = 'Present';
        end
        
        groupxs = [];
        for ci = 1:length(otherSideCorrToPlot)   
            otherCorrI =  find(valsByIndex.otherCorrects==otherSideCorrToPlot(ci));
            if otherSideCorrToPlot(ci)==0
                otherCorrLabs{ci} = 'Other incorr.';
            elseif otherSidePresToPlot(pi)==1
                otherCorrLabs{ci} = 'Other corr.';
            end
            
            bi = bi +1;
            if ci ~= 1 || (ci==1 && pi==1)
                bctr = bctr+bwid+bsep;
            else %add extra separation
                bctr = bctr+bwid+bsep*2.5;
            end
            
            bx=bctr+[-bwid/2 bwid/2];
            by=[ylims(1) ds(otherPresI,otherCorrI)];
            vertx=[bx; bx];
            verty=[by fliplr(by)];
            
            chandles(pi,ci)=fill(vertx(:), verty(:),otherCorrCols(ci,:),'EdgeColor',otherCorrCols(ci,:),'LineWidth',datLineWidth);
                        
            ctrs(bi)=bctr;
            groupxs = [groupxs bctr];
        end
        groupctrs = [groupctrs mean(groupxs)];
        

    end

    
    set(gca,'XTick',groupctrs,'XTickLabel',otherPresLabs,'LineWidth',axLineWidth);
    xlabel('Other side target');
    if si==1
        set(gca,'YTick',yticks,'YTickLabel',ytickLabs);
        ylabel('Ag');
    else
        set(gca,'YTick',yticks,'YTickLabel',{});
    end
    ylim(ylims);
    xlims = [0 ctrs(end)+bwid/2+bsep];
    
    if si == 1, legend(chandles(1,:),otherCorrLabs,'Location','NorthEast'); legend boxoff; end
    
    xlim(xlims);
    title(sideLabels{si});
end

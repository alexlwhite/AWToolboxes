function plotCuexSidexCongr(Ag,AgValsByIndex)

%Ag comes in with this structure:
%dim1 = single-task left; single-task right; dual-task
%dim2 = left stim; right stim; 
%dim3 = all; incongruent, congruent 

%re-organize into this form, to get rid of uncued sides 
%dim1 = single-task; dual-task
%dim2 = left; right
%dim3 = congruent, incongruent 

As = NaN(2,2,2); 

congs = 2:3;
for congi=1:2
    cong = congs(congi);
    for ci=1:3
        if ci<3
            As(1,ci,congi) = Ag(ci,ci,cong);
        else
            As(2,:,congi) = squeeze(Ag(ci,:,cong));
        end
    end
end

cueLabels = {'Single-task','Dual-task'}; 
sideLabels = AgValsByIndex.side;
congLabels = AgValsByIndex.congruency(congs);

%Colors
xColors = [46 63 153; 147 26 29]/255;
xHSV = rgb2hsv(xColors);

ylims = [0.5 1];
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

%Colors
cueXSideEdgeColors = zeros(2,2,2,3);
cueXSideFillColors = zeros(2,2,2,3);
for task = 1:1
    for cue = 1:2
        for cong = 1:2
            if cong==2 %CONGRUENT
                edgeHSV = xHSV(task,:).*[1 1 1.1]; %edge is usual color, 10% brigher
                
                if cue==2 %Focused
                    fillHSV = edgeHSV; %color fill
                else %divided: white fill
                    fillHSV = [0 0 1];
                end
            else %INCONGRUENT
                edgeHSV = xHSV(task,:).*[1 1 0.45]; %edge is slighlty dark than congruent
                if cue==2 %Focused: fill color
                    fillHSV =edgeHSV;
                else %divided:
                    %fillHSV = xHSV(x,:).*[1 0 0.96]; %gray fill
                    fillHSV = [0 0 1]; %white fill
                end
            end
            cueXSideEdgeColors(task,cue,cong,:) = hsv2rgb(edgeHSV);
            cueXSideFillColors(task,cue,cong,:) = hsv2rgb(fillHSV);
        end
    end
end

bwid=.7;
bsep=.3;

figure; hold on;
task = 1;
for ci = 1:2
    
    ds = squeeze(As(ci,:,:));
        
    subplot(1,2,ci);
    hold on;
    
    bi = 0;
    bctr = -bsep;
    ctrs = [];
    groupctrs = [];
    
    for si = 1:2
        groupxs = [];
        for oi = 1:2            
            bi = bi +1;
            if oi ~= 1 || (oi==1 && si==1)
                bctr = bctr+bwid+bsep;
            else %add extra separation
                bctr = bctr+bwid+bsep*2.5;
            end
            
            bx=bctr+[-bwid/2 bwid/2];
            by=[0 ds(si,oi)];
            vertx=[bx; bx];
            verty=[by fliplr(by)];
            
            chandles(si,oi)=fill(vertx(:), verty(:),squeeze(cueXSideFillColors(task,si,oi,:))','EdgeColor',squeeze(cueXSideEdgeColors(task,si,oi,:))','LineWidth',datLineWidth);
                        
            ctrs(bi)=bctr;
            groupxs = [groupxs bctr];
        end
        groupctrs = [groupctrs mean(groupxs)];
        

    end

    
    set(gca,'XTick',groupctrs,'XTickLabel',sideLabels,'LineWidth',axLineWidth);
    if ci==1
        set(gca,'YTick',yticks,'YTickLabel',ytickLabs);
        ylabel('Ag');
    else
        set(gca,'YTick',yticks,'YTickLabel',{});
    end
    ylim(ylims);
    xlims = [0 ctrs(end)+bwid/2+bsep];
    
    if ci == 1, legend(chandles(1,:),congLabels,'Location','NorthEast'); legend boxoff; end
    
    xlim(xlims);
    title(cueLabels{ci});
end
%set(gcf,'color','w','units','centimeters','pos',twoCondPlotSize);
%exportfig(gcf,fullfile(figFolder,sprintf('CuexSidexCong_M1_Task%i.eps',task)),'Format','eps','bounds','loose','color','rgb','LockAxes',0,'FontMode','fixed','FontSize',fSize);
end
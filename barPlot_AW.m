%% function barPlot_AW(ds, eb, opt)
% Alex White's function to make a bar plot. 
% 
% Inputs: 
% - ds: data values, in a 1 or 2D matrix. 
%       The first dimension defines groups that are separted more widely 
%       The second dimension defines values within a group 
% 
% - eb: Values for error bars on the bars. Same dimensions of ds. Set to
%       empty vector [] is none. 
% 
% - opt: structure with various plotting options. 
%       
function barPlot_AW(ds, eb, opt)

if ~isfield(opt,'barWidth')
    opt.barWidth = 0.2;
end
if ~isfield(opt,'level1Sep')
    opt.level1Sep = 0.6;
end
if ~isfield(opt,'level2Sep')
    opt.level2Sep = 0.25;
end
if ~isfield(opt,'xAxisMargin')
    opt.xAxisMargin = 0.3;
end
if ~isfield(opt,'edgeLineWidth')
    opt.edgeLineWidth = 1;
end
if ~isfield(opt,'errorBarWidth')
    opt.errorBarWidth = 1;
end

if ~isfield(opt,'ylims')
    drng=[min(ds(:)-eb(:)) max(ds(:)+eb(:))]; 
    opt.ylims = drng + [-1 1]*0.1*diff(drng);
end

if ~isfield(opt,'yticks')
    opt.yticks = linspace(opt.ylims(1), opt.ylims(2), 5);
end

if ~isfield(opt,'legendLoc')
    opt.legendLoc = 'NorthWest';
end


n1 = size(ds,1);
n2 = size(ds,2);

if ~isfield(opt,'fillColors')
    l2FillColors = hsv2rgb([linspace(0.2, 0.8, n2)' ones(n2,1)*0.8 ones(n2,1)*0.7]);
    opt.fillColors = zeros(n1,n2,3);
    for i1=1:n1
        opt.fillColors(i1,:,:) = l2FillColors;
    end
end

if ~isfield(opt,'edgeColors')
    opt.edgeColors = opt.fillColors;
end

if ~isfield(opt,'errorBarColors')
    opt.errorBarColors = opt.edgeColors;
end

%Set bar centers
ctr = 0;
barCenters = zeros(n1,n2);

for i1 = 1:n1
    for i2 = 1:n2
        if i1==1 && i2==1
            ctr = ctr + opt.xAxisMargin;
        elseif i2==1
            ctr = ctr + opt.level1Sep;
        else
            ctr = ctr + opt.level2Sep;
        end
        barCenters(i1,i2) = ctr;
    end
end

xlims = [0 ctr+opt.xAxisMargin];

figure; hold on;
handles = zeros(1,n2);

if prod(opt.ylims)<0
    plot(xlims,[0 0],'k-');
end

for i1 = 1:n1
    for i2 = 1:n2
        bx = barCenters(i1,i2)+[-0.5 0.5]*opt.barWidth;
        by = [0 ds(i1, i2)]; 
        vertx=[bx; bx];
        verty=[by fliplr(by)];
       
        %bar
        handles(i2) = fill(vertx(:), verty(:),squeeze(opt.fillColors(i1,i2,:))','EdgeColor',squeeze(opt.edgeColors(i1,i2,:))','LineWidth',opt.edgeLineWidth);
        
        %error bar
        if ~isempty(eb)
            plot([1 1]*barCenters(i1,i2),ds(i1,i2)+[-1 1]*eb(i1,i2),'-','Color',squeeze(opt.errorBarColors(i1,i2,:)),'LineWidth',opt.errorBarWidth);
        end
    end
end

xlim(xlims); ylim(opt.ylims);

l1Centers = mean(barCenters,2); 
set(gca,'XTick',l1Centers,'YTick',opt.yticks); 
if isfield(opt,'xTickLabs')
    set(gca,'XTickLabel',opt.xTickLabs); 
end

if isfield(opt,'xLab')
    xlabel(opt.xLab);
end
if isfield(opt,'yLab')
    ylabel(opt.yLab); 
end

if isfield(opt,'legendLabs')
    legend(handles,opt.legendLabs,'Location',opt.legendLoc);
end



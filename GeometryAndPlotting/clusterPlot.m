%% function  [groupCenters, opt] = clusterPlot_AW(ds, opt)
% Alex White's function to make a version of a bar plot that has individual data points in clusters. 
% 
% Inputs: 
% - ds: cell array of individual data points, with one cell per condition.  
%       The first dimension defines groups that are separted more widely 
%       The second dimension defines values within a group 
% 
% 
% - opt: structure with various plotting options: 
%    barWidth
%    meanLineWidth
%    errorBarWidth
%    fillColors
%    edgeColors
%    errorBarColors
%    level1Sep
%    level2Sep
%    xAxisMargin
%    xLab
%    xTickLabs
%    ylims
%    yLab
%    doYLab
%    yticks
%    doLegend
%    legendLabs
%    legendLoc
%    legendTitle
%    lev1ForLegend (which value of level 1 to use for handles for legend)
% 
% Outputs: 
% - barCenters: NxM matrix of bar centers, where N is the size of dimension
% 1 in the data and M is the size of dimension 2
% 

function [barCenters, opt] = clusterPlot_AW(ds, opt)

n1 = size(ds,1);
n2 = size(ds,2);
%check data size: each element of ds should be a column vector 
for i1=1:n1
    for i2=1:n2
        ds{i1,i2} = squeeze(ds{i1,i2});
        dsz=size(ds{i1,i2});
        if length(dsz)>2 || all(dsz>1)
            error('Each cell of input ds must be a column vector')
        end
        if dsz(1)==1
            ds{i1,i2} = ds{i1,i2}';
        end
    end
end
        

if ~isfield(opt,'markSz')
    opt.markSz = 10;
end
if ~isfield(opt, 'symbols')
    opt.symbols = repmat({'o'}, n1, n2);
end
if ~isfield(opt, 'barWidth')
    opt.barWidth = 0.15;
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
if ~isfield(opt,'meanLineWidth')
    opt.meanLineWidth = 3;
end
if ~isfield(opt, 'nVertBands')
    opt.nVertBands = 10;
end
if ~isfield(opt, 'doErrorBar')
    opt.doErrorBar = true;
end
if ~isfield(opt, 'symmetricErrorBar')
    opt.symmetricErrorBar = true;
end

if ~isfield(opt,'errorBarWidth')
    opt.errorBarWidth = 1;
end
if ~isfield(opt, 'legendTitle')
    opt.legendTitle = '';
end
if ~isfield(opt,'errorBarColors')
    opt.errorBarColors = zeros(size(ds,1), size(ds,2), 3);
end

if ~isfield(opt, 'connectIndivPts')
    opt.connectIndivPts = true;
end
if ~isfield(opt, 'indivPtConnectColor')
    opt.indivPtConnectColor = 0.8*ones(1,3);
end

if ~isfield(opt,'doYLab')
    opt.doYLab = true;
end

if ~isfield(opt,'doYTickLab')
    opt.doYTickLab = true;
end
if ~isfield(opt,'doLegend')
    opt.doLegend = true;
end
if ~isfield(opt,'legendLoc')
    opt.legendLoc = 'NorthWest';
end
if ~isfield(opt,'lev1ForLegend')
    opt.lev1ForLegend = size(ds,1);
end


if ~isfield(opt,'fillColors')
    l2FillColors = hsv2rgb([linspace(0.3, 0.8, n2)' ones(n2,1)*0.8 ones(n2,1)*0.7]);
    opt.fillColors = zeros(n1,n2,3);
    for i1=1:n1
        opt.fillColors(i1,:,:) = l2FillColors;
    end
elseif n2==1 && size(opt.fillColors,3)==1 %reshape fillColors 
    opt.fillColors = reshape(opt.fillColors,[n1 n2 3]);
end

if ~isfield(opt,'edgeColors')
    opt.edgeColors = opt.fillColors;
elseif n2==1 && size(opt.edgeColors,3)==1 %reshape edgeColors 
    opt.edgeColors = reshape(opt.edgeColors,[n1 n2 3]);
end

if ~isfield(opt, 'meanLineColors')
    opt.meanLineColors = opt.edgeColors*0.8;
end

%Set bar centers and collect all points 
ctr = 0;
barCenters = zeros(n1,n2);
allPts = [];
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
        allPts = [allPts; ds{i1,i2}(:)];
    end
end

xlims = [0 ctr+opt.xAxisMargin];

yrng = [min(allPts) max(allPts)];
if ~isfield(opt,'ylims')
    opt.ylims = yrng + [-1 1]*0.1*diff(yrng);
elseif opt.ylims(1)>yrng(1) || opt.ylims(2)<yrng(2)
    fprintf(1,'\n(%s) Input ylims excluded some data. Resetting to not do that\n', mfilename);
    opt.ylims = yrng + [-1 1]*0.1*diff(yrng);
end

if ~isfield(opt,'yticks')
    nTicks = 5;
    tickDiff = diff(opt.ylims)/(nTicks-1);
    
    orderOfMag = round(log10(tickDiff));
    
    logDiff = 1 - orderOfMag;
    ts = tickDiff*10^logDiff;
    ts = round(ts);
    tickDiff = ts/(10^logDiff);
    
    startTick = tickDiff*round(opt.ylims(1)/tickDiff);
    
    opt.yticks = startTick:tickDiff:opt.ylims(2);
end    



%divide data points into vertical "bands". Dots of the same condition in
%the same band get shifted horizontally relative to each other 
[~, edges] = histcounts(allPts, opt.nVertBands);
%nudge the last egde up to avoid errors later 
edges(end)=edges(end)*1.001; 
%figure out the max number of points in any band across individual
%conditions
maxN = 0;
for i1 = 1:n1
    for i2 = 1:n2
        thisN = histcounts(ds{i1,i2}, edges);
        maxN = max([maxN max(thisN)]);
    end
end

pointHorizSep = 0.9*(opt.barWidth)/(maxN);

%assign each dot an x and y coordinate 
allX = cell(size(ds));
for i1 = 1:n1
    for i2 = 1:n2
        allX{i1,i2}=NaN(size(ds{i1,i2}));
        thisN = histcounts(ds{i1,i2}, edges);
        for ei=1:opt.nVertBands
            if thisN(ei)>0
                xis = edges(ei)<=ds{i1,i2} & ds{i1,i2}<edges(ei+1);
                wid = thisN(ei)*pointHorizSep;
                try
                    allX{i1,i2}(xis) = barCenters(i1,i2)+linspace(-wid/2, wid/2, thisN(ei));
                catch
                    keyboard
                end
            end
        end
    end
end


hold on;
handles = zeros(n1,n2);

if prod(opt.ylims)<0
    plot(xlims,[0 0],'k-');
end

%if requesed, plot lines connecting matched data points 
if opt.connectIndivPts 
   ys = []; xs = [];
   for i1=1:n1
       ys = [ys cell2mat(ds(i1,:))];
       %xs = [xs barCenters(i1,:)];
       xs = [xs cell2mat(allX(i1,:))];
   end
   
   for pti=1:size(ys,1)
       plot(xs(pti,:), ys(pti,:), '-','Color',opt.indivPtConnectColor);
   end
end

for i1 = 1:n1
    for i2 = 1:n2
        bx = barCenters(i1,i2)+[-0.5 0.5]*opt.barWidth;
        by = ones(1,2)*nanmean(ds{i1, i2}); 
      
       
        %dots
        handles(i1,i2)  = plot(allX{i1,i2}, ds{i1,i2}, opt.symbols{i1,i2}, 'MarkerFaceColor',squeeze(opt.fillColors(i1,i2,:))', 'MarkerEdgeColor', squeeze(opt.edgeColors(i1,i2,:))');
        
        %line for mean
        plot(bx, by, '-', 'Color', squeeze(opt.meanLineColors(i1,i2,:)), 'LineWidth',opt.meanLineWidth);
                
        %error bar
        if opt.doErrorBar
            if opt.symmetricErrorBar
                plot([1 1]*barCenters(i1,i2), nanmean(ds{i1,i2})+[-1 1]*standardError(ds{i1,i2}), '-', 'Color', squeeze(opt.errorBarColors(i1,i2,:)), 'LineWidth', opt.errorBarWidth);
            else
                plot([1 1]*barCenters(i1,i2), boyntonBootstrap(@mean,ds{i1,i2},1000,68.27,true),'-','Color', squeeze(opt.errorBarColors(i1,i2,:)), 'LineWidth', opt.errorBarWidth);
            end
        end
    end
end


l1Centers = mean(barCenters,2); 
set(gca,'XTick',l1Centers,'YTick',opt.yticks); 
if isfield(opt,'xTickLabs')
    set(gca,'XTickLabel',opt.xTickLabs); 
end
if ~opt.doYTickLab
    set(gca,'YTickLabel',{});
end
if isfield(opt,'xLab')
    xlabel(opt.xLab);
end
if isfield(opt,'yLab') && opt.doYLab
    ylabel(opt.yLab); 
end

xlim(xlims); 
ylim(opt.ylims);

if isfield(opt,'legendLabs') && opt.doLegend
    leg = legend(handles(opt.lev1ForLegend,:),opt.legendLabs,'Location',opt.legendLoc,'AutoUpdate','off');

    if ~isempty(opt.legendTitle)
        title(leg, opt.legendTitle);
    end
end



%% function  [barCenters, opt, legendIcons] = clusterPlot_AW(ds, opt)
% Alex White's function to make a version of a bar plot that has individual data points in clusters.
%
% Inputs:
% - ds: cell array of individual data points, with one cell per condition.
%       The first dimension defines groups that are separted more widely
%       The second dimension defines values within a group
%
%
% - opt: structure with various plotting options: NEEDS UPDATING
%    fillColors
%    edgeColors
%    barWidth
%    meanLineWidth
%    connectLev1IndivPts
%    connectLev2IndivPts
%    indivPtConnectColor
%    errorBarWidth
%    doErrorBar
%    errorBarColors
%    errorBarCI
%    errorBarType
%    symmetricErrorBar
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
%    markSz 
%    meanSymbol
%    meanLineWidth
%    meanDotSize
%    nVertBands
% Outputs:
% - barCenters: NxM matrix of bar centers, where N is the size of dimension
% 1 in the data and M is the size of dimension 2
%

function [barCenters, opt, legendIcons] = clusterPlot(ds, opt)

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
    opt.markSz = 100;
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
if ~isfield(opt, 'meanSymbol')
    opt.meanSymbol = 'd'; %'-' means a line
end
if ~isfield(opt,'meanLineWidth')
    opt.meanLineWidth = 4;
end
if ~isfield(opt,'meanDotSize')
    opt.meanDotSize = 10;
end

if ~isfield(opt, 'nVertBands')
    opt.nVertBands = 50;
end

if ~isfield(opt, 'doErrorBar')
    opt.doErrorBar = true;
end
if ~isfield(opt, 'symmetricErrorBar')
    opt.symmetricErrorBar = true;
end
if ~isfield(opt,'errorBarCI')
    opt.errorBarCI = 68.27;
end

if ~isfield(opt, 'errorBarType')
    %possibilities: 'box', 'line'
    opt.errorBarType = 'line';
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

if ~isfield(opt, 'connectLev1IndivPts')
    opt.connectLev1IndivPts = false;
end
if ~isfield(opt, 'connectLev2IndivPts')
    opt.connectLev2IndivPts = true;
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

if ~isfield(opt, 'meanColors')
    opt.meanColors = opt.edgeColors*0.7;
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
    nTicks = 6;
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

pointHorizSep = 0.5*(opt.barWidth)/(maxN);

%assign each dot an x and y coordinate
allX = cell(size(ds));
for i1 = 1:n1
    if opt.connectLev2IndivPts %save info about x-pos and y-bands so we can make sure each individual's points are at the same x-position
        connectXJitts = zeros(size(ds{1,1}));
        connectYBands = zeros(size(connectXJitts));
    end
    for i2 = 1:n2
        allX{i1,i2}=NaN(size(ds{i1,i2}));
        thisN = histcounts(ds{i1,i2}, edges);
        for ei=1:opt.nVertBands
            if thisN(ei)>0
                xis = edges(ei)<=ds{i1,i2} & ds{i1,i2}<edges(ei+1);
                wid = thisN(ei)*pointHorizSep;
                if thisN(ei)==1
                    allX{i1,i2}(xis) = barCenters(i1,i2);
                else
                    allX{i1,i2}(xis) = barCenters(i1,i2)+linspace(-wid/2, wid/2, thisN(ei));
                end
                if opt.connectLev2IndivPts
                    connectYBands(xis, i2) = ei;
                    connectXJitts(xis, i2) = allX{i1,i2}(xis)-barCenters(i1,i2);
                end
            end
        end
    end
    
    % make sure that all the pairs of connected points are at the same
    % relative x-position
    if opt.connectLev2IndivPts && n2==2 %currently only works when there are 2 levels of dim2, that is, just connecting pairs of points. See use of "dim2adjust"
        maxRep = 10000;
        repI = 0;
        allGood = false;
        dimToAdjust = 1;
        while ~allGood && repI<maxRep
            repI = repI+1;
            alreadyFixed = [];
            dimToAdjust = 3-dimToAdjust;
            ranAground = false;
            xOrder = randperm(size(connectXJitts,1));
            for ii = xOrder
                needsFix = std(connectXJitts(ii,:), 0, 2)>10^-10;
                if needsFix && ~ranAground
                    
                    badXJitt = connectXJitts(ii, dimToAdjust);
                    goalXJitt = connectXJitts(ii, 3-dimToAdjust);
                    
                    thisYBand = connectYBands(ii,dimToAdjust);
                    inThisYBand = find(connectYBands(:,dimToAdjust)==thisYBand);
                    theseYIs = setdiff(inThisYBand, alreadyFixed);
                    if length(theseYIs)<2
                        if length(inThisYBand)==1  
                            connectXJitts(ii, dimToAdjust) = goalXJitt;
                            alreadyFixed = [alreadyFixed ii];
                        else
                            ranAground = true;
                        end
                    else
                        theseXJitts = connectXJitts(theseYIs, dimToAdjust);
                        inTheWay = find(theseXJitts==goalXJitt);
                        thisII = find(theseYIs==ii);
                        
                        if ~isempty(inTheWay)
                            %swap the problematic x-jitter with the other one that has the needed value:
                            
                            newOrder = 1:length(theseYIs);
                            swapI = inTheWay(1);
                            newOrder(swapI) = length(newOrder)+1;
                            newOrder(thisII) = swapI;
                            newOrder(swapI) = thisII;
                            newXJitts = theseXJitts(newOrder);
                            %change the swapped one to not fixed
                            alreadyFixed = setdiff(alreadyFixed, theseYIs(swapI));
                            
                        else %if none in this y-band have the required x-jitt, just set it to that
                            newXJitts = theseXJitts;
                            newXJitts(thisII) = goalXJitt;
                        end
                        
                        connectXJitts(theseYIs, dimToAdjust) = newXJitts;
                        alreadyFixed = [alreadyFixed ii];
                        gotFixed = std(connectXJitts(ii,:), 0, 2)==0;
                        if ~gotFixed
                            keyboard
                        end
                    end
                    
                end
            end
            %check
            allGood = all(abs(diff(connectXJitts,1,2))<10^-15);
        end
        
        
        %re-set the original x-positions
        for i2 = 1:n2
            allX{i1,i2} = barCenters(i1,i2)+connectXJitts(:,i2);
        end
        %double cheeck
        theseX = cell2mat(allX(i1,:));
        itWorked = std(diff(theseX,1,2))<10^-15;
        if ~itWorked
            keyboard
        end
    end
end


hold on;
handles = zeros(n1,n2);

if prod(opt.ylims)<0
    plot(xlims,[0 0],'k-');
end

%if requested, plot lines connecting matched data points
if opt.connectLev2IndivPts
    for i1=1:n1
        xs = cell2mat(allX(i1,:));
        ys = cell2mat(ds(i1,:));
        for pti=1:size(ys,1)
            plot(xs(pti,:), ys(pti,:), '-','Color',opt.indivPtConnectColor);
        end
    end
end

if opt.connectLev1IndivPts
%     ys = []; xs = [];
%     for i1=1:n1
%         ys = [ys cell2mat(ds(i1,:))];
%         xs = [xs cell2mat(allX(i1,:))];
%     end
%     
%     for pti=1:size(ys,1)
%         plot(xs(pti,:), ys(pti,:), '-','Color',opt.indivPtConnectColor);
%     end

    for i1=1:(n1-1)
        x1 = squeeze(allX{i1, n2, :});
        x2 = squeeze(allX{i1+1, 1, :});
        
        xs = [x1 x2];
        
        y1 = squeeze(ds{i1, n2, :});
        y2 = squeeze(ds{i1+1, 1, :});
        
        ys = [y1 y2];
        for pti=1:size(ys,1)
            plot(xs(pti,:), ys(pti,:), '-','Color',opt.indivPtConnectColor);
        end
        
        
    end

    
    fprintf(1,'\n\n(%s) WARNNING: connectLev1IndivPts is not coded to ensure connected points have same-xjitter so slopes are comparable\n', mfilename);
end


for i1 = 1:n1
    for i2 = 1:n2
        
        
        %dots
        hdot  = scatter(allX{i1,i2}, ds{i1,i2});
        hdot.Marker = opt.symbols{i1,i2};
        hdot.SizeData = opt.markSz;
        hdot.MarkerFaceColor = squeeze(opt.fillColors(i1,i2,:))';
        hdot.MarkerEdgeColor = squeeze(opt.edgeColors(i1,i2,:))';
        
        if isfield(opt, 'dotFaceAlpha')
            hdot.MarkerFaceAlpha = opt.dotFaceAlpha;
        end
        
        handles(i1,i2) = hdot;
        
        if strcmp(opt.meanSymbol, '-')
            %line for mean
            bx = barCenters(i1,i2)+[-0.5 0.5]*opt.barWidth;
            by = ones(1,2)*nanmean(ds{i1, i2});
            plot(bx, by, '-', 'Color', squeeze(opt.meanColors(i1,i2,:)), 'LineWidth',opt.meanLineWidth);
        else
            plot(barCenters(i1,i2), nanmean(ds{i1, i2}), opt.meanSymbol, 'Color', squeeze(opt.meanColors(i1,i2,:)), 'MarkerSize',opt.meanDotSize);
        end
        %error bar as a box!
        if opt.doErrorBar
            if opt.symmetricErrorBar
                if strcmp(opt.errorBarType, 'line')
                    plot([1 1]*barCenters(i1,i2), nanmean(ds{i1,i2})+[-1 1]*standardError(ds{i1,i2}), '-', 'Color', squeeze(opt.errorBarColors(i1,i2,:)), 'LineWidth', opt.errorBarWidth);
                elseif strcmp(opt.errorBarType, 'box')
                    rectWid = opt.barWidth;
                    rectHei = 2*standardError(ds{i1,i2});
                    rectX = barCenters(i1,i2)-rectWid/2; %lower left corner
                    rectY = nanmean(ds{i1,i2}-rectHei/2);
                    
                    rectangle('Position', [rectX rectY rectWid rectHei],'EdgeColor', [0 0 0]);
                end
            else
                theseDats = ds{i1,i2};
                theseDats = theseDats(~isnan(theseDats));
                CI = boyntonBootstrap(@mean,theseDats,1000,opt.errorBarCI,true);
                
                if strcmp(opt.errorBarType, 'line')
                    plot([1 1]*barCenters(i1,i2), CI,'-','Color', squeeze(opt.errorBarColors(i1,i2,:)), 'LineWidth', opt.errorBarWidth);
                elseif strcmp(opt.errorBarType, 'box')
                    rectWid = opt.barWidth;
                    rectX = barCenters(i1,i2)-rectWid/2;
                    rectHei = abs(diff(CI));
                    rectY = min(CI);
                    rectangle('Position', [rectX rectY rectWid rectHei],'EdgeColor', [0 0 0]);
                end
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
    [leg, legendIcons] = legend(handles(opt.lev1ForLegend,:),opt.legendLabs,'Location',opt.legendLoc,'AutoUpdate','off');
    
    if ~isempty(opt.legendTitle)
        title(leg, opt.legendTitle);
    end
else
    legendIcons = [];
end



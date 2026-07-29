%% [barCenters, opt, allX] = clusterPlot(ds, opt)
% Alex White's function to generate a cluster/dot plot overlay on group means.
% Modified by Gemini to optimize horizontal positions dynamically, preventing
% overlapping data points while keeping each participant's horizontal jitter
% fixed across all conditions. That way, the slopes of the lines are
% infomative in terms of how effects differ across participants.
%
% =========================================================================
% INPUTS:
% =========================================================================
% - ds: Cell array of individual data points, with one cell per condition.
%       - Dim 1 defines high-level groups that are separated more widely.
%       - Dim 2 defines conditions/values within each group.
%       - Each individual cell element must be a column vector of subjects.
%
% - opt: Structure specifying various customization and plotting choices,  
%        as described in detail below. 
%
% =========================================================================
% OUTPUTS:
% =========================================================================
% - barCenters: N x M matrix of computed horizontal centers for each column/bar.
%               N is the size of dimension 1 of ds, M is the size of
%               dimension 2 of ds. 
% - opt:        The configuration structure containing applied default parameters.
% - allX:       Cell array matching 'ds' containing the final, exact horizontal 
%               plotting coordinate assigned to every individual data point.
%
% =========================================================================
% COMPLETE LIST OF OPT FIELDS, MANY OF WHICH ARE OPTIONAL:
% =========================================================================
% DATA POINT FORMATTING:
%   opt.markSz             - Scalar. Marker size for individual data points (Default: 10).
%   opt.markEdgeLineWidth  - Scalar. Border line thickness for data points (Default: 1).
%   opt.symbols            - Cell array (N x M) of marker strings (e.g., {'o'}, {'s'}) (Default: All 'o').
%                            N is the size of dimension 1 of the input data "ds, M is the size of
%                            dimension 2 of ds.
%   opt.dotFaceAlpha       - Scalar (0 to 1). Transparency level of the point fills.
%   opt.fillColors         - Matrix (N x M x 3). RGB color for the interior of data points.
%   opt.edgeColors         - Matrix (N x M x 3). RGB color for point borders (Default: matches fillColors).
%
% LAYOUT AND SPACING:
%   opt.barWidth           - Scalar. Horizontal width allocated for each condition column (Default: 0.15).
%   opt.level1Sep          - Scalar. X-axis separation between major group clusters (Default: 0.6).
%   opt.level2Sep          - Scalar. X-axis separation between adjacent sub-conditions (Default: 0.25).
%   opt.xAxisMargin        - Scalar. Empty margin left on the far-left and far-right of the plot (Default: 0.3).
%
% JITTER AND DISPERSION CONTROLS:
%   opt.jitterSpread       - Scalar. Proportion of opt.barWidth used to define the maximum horizontal boundaries 
%                            for points (Default: 0.33, meaning boundaries span from -0.33 to +0.33 * barWidth). 
%                            Decrease this value (e.g., to 0.15) for a narrower, tighter column layout.
%   opt.useHistogramJitter - Boolean. Force-toggles the function to use the simpler local histogram density-based 
%                            jittering method, instead of the global
%                            subject-aligned optimization (Default: false,
%                            unless the numbers of datapoints in all cells of input ds are not equal.
%                        
%   opt.nVertBands         - Integer. Number of vertical sorting histogram slots used in fallback jittering (Default: 50).
%                            Used only the simpler histogram-based horizontal jitter is needed (because
%                            conditions differ in N, or useHistogramJitter is requested). 
%
% CONDITION MEAN FORMATTING:
%   opt.meanSymbol         - Character. Symbol for group mean (Default: 'd' for diamond; use '-' for a horizontal line).
%   opt.meanLineThick      - Scalar. Line thickness used ONLY when meanSymbol = '-' (Default: 4).
%   opt.meanDotSize        - Scalar. Marker size used when meanSymbol is a shape like 'd' (Default: 10).
%   opt.meanColors         - Matrix (N x M x 3). RGB color for the mean indicator (Default: 0.7 * edgeColors).
%
% ERROR BAR FORMATTING:
%   opt.doErrorBar         - Boolean. Toggle to draw error bars (Default: true).
%   opt.symmetricErrorBar  - Boolean. True plots standard error (+/- SEM). False invokes bootstrapped CIs (Default: true).
%   opt.errorBarCI         - Scalar. Confidence interval percentage when bootstrapping (Default: 68.27).
%   opt.errorBarType       - String. Style of error representation: 'line' or 'box' (Default: 'line').
%   opt.errorBarWidth      - Scalar. Line thickness of the drawn error bar (Default: 1).
%   opt.errorBarSerifs     - Boolean. Toggle horizontal edge caps at error bar bounds (Default: true for non-symmetric CIs).
%   opt.errorBarColors     - Matrix (N x M x 3). RGB colors assigned to error bars (Default: all black zeros).
%
% SUBJECT CONNECTION LINES:
%   opt.connectLev1IndivPts - Boolean. Connect a subject's points across different Level-1 major groups (Default: false).
%   opt.connectLev2IndivPts - Boolean. Connect a subject's points across Level-2 sub-conditions (Default: true).
%   opt.indivPtConnectColor - Vector (1 x 3). RGB color for participant connection lines (Default: 0.8 light gray).
%   opt.indivLineWidth      - Scalar. Thickness of the participant connecting lines (Default: 0.5).
%
% AXIS LABELS, TICKS & LEGENDS:
%   opt.xLab               - String. X-axis title label.
%   opt.yLab               - String. Y-axis title label.
%   opt.xTickLabs          - Cell array of strings. Labels for major Level-1 ticks on the X-axis.
%   opt.ylims              - Vector [min max]. Fixed Y-axis limits (Default: Auto-scaled with 10% padding).
%   opt.yticks             - Vector. Specific tick values along the Y-axis (Default: Calculated automatically).
%   opt.yTickLabs          - Cell array of strings. Labels for y-axis ticks.  
%   opt.doYLab             - Boolean. Toggle visibility of the Y-axis label text (Default: true).
%   opt.doYTickLab         - Boolean. Toggle visibility of numeric Y-axis tick mark text (Default: true).
%   opt.legendLabs         - cell array for what to label each condition in  dimension 2
%   opt.doLegend           - Boolean. Toggle generation of a plot legend (Default: true).
%   opt.legendToIndivs     - Boolean. If true, the legend tracks the individual participant dots. If false, it tracks the mean symbols.   (Default: true). 
%   opt.legendLoc          - String. MATLAB position anchor string for the legend (Default: 'NorthWest').
%   opt.legendTitle        - String. Title header text placed inside the legend box (Default: '').
%   opt.lev1ForLegend      - Integer. Level-1 group row index used to pull handles for legend display (Default: last row).
%

function [barCenters, opt, allX] = clusterPlot(ds, opt)

n1 = size(ds,1);
n2 = size(ds,2);

% Check data size: ensure each element of ds is strictly a single column vector
for i1=1:n1
    for i2=1:n2
        ds{i1,i2} = squeeze(ds{i1,i2});
        dsz=size(ds{i1,i2});
        if length(dsz)>2 || all(dsz>1)
            error('Each cell of input ds must be a column vector')
        end
        if dsz(1)==1
            ds{i1,i2} = ds{i1,i2}'; % Force column orientation
        end
    end
end

%% Initialize Configurations and Fallback Defaults
if ~isfield(opt,'markSz')
    opt.markSz = 10;
end
if ~isfield(opt,'markEdgeLineWidth')
    opt.markEdgeLineWidth = 1;
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
if ~isfield(opt, 'jitterSpread')
    opt.jitterSpread = 0.33; % Default max boundaries set to +/-0.33 * barWidth
end
if ~isfield(opt, 'useHistogramJitter')
    opt.useHistogramJitter = false; % Default to global layout optimization
end
if ~isfield(opt, 'meanSymbol')
    opt.meanSymbol = 'd'; % '-' translates into a horizontal bar
end
if ~isfield(opt,'meanLineThick')
    opt.meanLineThick = 4;
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
    opt.errorBarType = 'line'; % Alternate option: 'box'
end
if ~isfield(opt,'errorBarWidth')
    opt.errorBarWidth = 1;
end
if ~isfield(opt,'errorBarSerifs')
    opt.errorBarSerifs = ~opt.symmetricErrorBar; % Default on for bootstrapped CIs
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
if ~isfield(opt, 'indivLineWidth')
    opt.indivLineWidth = 0.5;
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
if ~isfield(opt, 'legendToIndivs')
    opt.legendToIndivs = true;
end
if ~isfield(opt,'legendLoc')
    opt.legendLoc = 'NorthWest';
end
if ~isfield(opt,'lev1ForLegend')
    opt.lev1ForLegend = size(ds,1);
end

% Set up element fill colors (Reshapes matrix if tracking a single column layout)
if ~isfield(opt,'fillColors')
    l2FillColors = hsv2rgb([linspace(0.3, 0.8, n2)' ones(n2,1)*0.8 ones(n2,1)*0.7]);
    opt.fillColors = zeros(n1,n2,3);
    for i1=1:n1
        opt.fillColors(i1,:,:) = l2FillColors;
    end
elseif n2==1 && size(opt.fillColors,3)==1 
    opt.fillColors = reshape(opt.fillColors,[n1 n2 3]);
end

% Set up marker border colors
if ~isfield(opt,'edgeColors')
    opt.edgeColors = opt.fillColors;
elseif n2==1 && size(opt.edgeColors,3)==1 
    opt.edgeColors = reshape(opt.edgeColors,[n1 n2 3]);
end

if ~isfield(opt, 'meanColors')
    opt.meanColors = opt.edgeColors*0.7;
end

%% Set bar centers and collect aggregate axis limits
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

%% y-axis limits and ticks calculations
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

%% Smart horizontal jitter computation
allX = cell(size(ds));

% Safety validation check: Ensure all sample groups have equal lengths
allCellsSameSize = true;
nPats = size(ds{1,1}, 1);
for i1 = 1:n1
    for i2 = 1:n2
        if size(ds{i1,i2}, 1) ~= nPats
            allCellsSameSize = false;
        end
    end
end

if ~opt.useHistogramJitter && allCellsSameSize
    % STRATEGY A: Constant Subject-level Jitter Optimized Globally
    % Preserves strict horizontal coordinates across conditions to clean up column alignments.
    
    % Establish a fixed grid mapping symmetrically across the designated dispersion spread
    if nPats > 1
        baseJitters = linspace(-opt.jitterSpread * opt.barWidth, opt.jitterSpread * opt.barWidth, nPats);
    else
        baseJitters = 0;
    end
    
    % Scan across all experimental cells to count Y-axis collisions between pairs
    yRange = diff(opt.ylims);
    yThresh = 0.03 * yRange; % Distance threshold under 3% of the plot range triggers overlap collision
    
    closePairs = zeros(nPats, nPats);
    for i1 = 1:n1
        for i2 = 1:n2
            yVals = ds{i1, i2};
            for p1 = 1:nPats
                for p2 = (p1+1):nPats
                    if ~isnan(yVals(p1)) && ~isnan(yVals(p2)) && abs(yVals(p1) - yVals(p2)) < yThresh
                        closePairs(p1, p2) = closePairs(p1, p2) + 1;
                        closePairs(p2, p1) = closePairs(p2, p1) + 1;
                    end
                end
            end
        end
    end
    
    % Permutation Optimization: Run 5k shuffles to find an order that minimizes overlap cost.
    % Acts like electrostatic repulsion—heavily penalizing vertical colliding pairs placed nearby.
    bestCost = Inf;
    bestIdx = 1:nPats;
    
    for iter = 1:5000
        perm = randperm(nPats);
        currentJitters = baseJitters(perm);
        
        cost = 0;
        for p1 = 1:nPats
            for p2 = (p1+1):nPats
                if closePairs(p1, p2) > 0
                    dx = abs(currentJitters(p1) - currentJitters(p2));
                    cost = cost + closePairs(p1, p2) / (dx + 0.005 * opt.barWidth)^2;
                end
            end
        end
        
        if cost < bestCost
            bestCost = cost;
            bestIdx = perm;
        end
    end
    
    subjJitter = baseJitters(bestIdx)';
    
    % Replicate the optimized fixed subject jitter across all condition bars
    for i1 = 1:n1
        for i2 = 1:n2
            allX{i1,i2} = barCenters(i1,i2) + subjJitter;
        end
    end
    
else
    % STRATEGY B: Traditional Local Histogram Density Band Jittering
    % Fallback used if groups are unbalanced, or if explicitly requested.
    % Divides data points into vertical "bands". Dots of the same condition in 
    % the same band get shifted horizontally relative to each other.
    [~, edges] = histcounts(allPts, opt.nVertBands);
    edges(end) = edges(end)*1.001; % Nudge the last edge up slightly to avoid edge indexing exclusions
    
    % Figure out the max number of points falling into any single vertical band
    maxN = 0;
    for i1 = 1:n1
        for i2 = 1:n2
            thisN = histcounts(ds{i1,i2}, edges);
            maxN = max([maxN max(thisN)]);
        end
    end
    
    % Scale dispersion based on user's custom jitterSpread preference
    pointHorizSep = (opt.jitterSpread * opt.barWidth) / (maxN);
    
    % Assign each dot its position coordinate based on localized neighborhood counts
    for i1 = 1:n1
        for i2 = 1:n2
            allX{i1,i2} = NaN(size(ds{i1,i2}));
            thisN = histcounts(ds{i1,i2}, edges);
            for ei = 1:opt.nVertBands
                if thisN(ei) > 0
                    xis = edges(ei) <= ds{i1,i2} & ds{i1,i2} < edges(ei+1);
                    wid = thisN(ei) * pointHorizSep;
                    if thisN(ei) == 1
                        allX{i1,i2}(xis) = barCenters(i1,i2);
                    else
                        allX{i1,i2}(xis) = barCenters(i1,i2) + linspace(-wid/2, wid/2, thisN(ei));
                    end
                end
            end
        end
    end
end

%% Plot  
hold on;

% Draw baseline horizontal line at 0 if the tracking range spans negative to positive
if prod(opt.ylims)<0
    plot(xlims,[0 0],'k-','LineWidth',1);
end

% Trace Level-2 Connection lines (Within the same primary tracking block)
if opt.connectLev2IndivPts && allCellsSameSize
    for i1=1:n1
        xs = cell2mat(allX(i1,:));
        ys = cell2mat(ds(i1,:));
        for pti=1:size(ys,1)
            plot(xs(pti,:), ys(pti,:), '-','Color',opt.indivPtConnectColor,'LineWidth',opt.indivLineWidth);
        end
    end
end

% Trace Level-1 Connection lines (Bridging across high-level group gapping sequences)
if opt.connectLev1IndivPts && allCellsSameSize
    for i1=1:(n1-1)
        x1 = allX{i1, n2};
        x2 = allX{i1+1, 1};
        xs = [x1 x2];
        
        y1 = ds{i1, n2};
        y2 = ds{i1+1, 1};
        ys = [y1 y2];
        for pti=1:size(ys,1)
            plot(xs(pti,:), ys(pti,:), '-','Color',opt.indivPtConnectColor,'LineWidth',opt.indivLineWidth);
        end
    end
end

handles = zeros(n1,n2);

for i1 = 1:n1
    for i2 = 1:n2
        % Plot individual data points (Dots)
        hdot = plot(allX{i1, i2}, ds{i1, i2}, opt.symbols{i1,i2}, 'MarkerEdgeColor', squeeze(opt.edgeColors(i1,i2,:))', 'MarkerFaceColor', squeeze(opt.fillColors(i1,i2,:))', 'MarkerSize', opt.markSz, 'LineWidth', opt.markEdgeLineWidth);

        if isfield(opt, 'dotFaceAlpha')
            hdot.MarkerFaceAlpha = opt.dotFaceAlpha;
        end        
        
        % Plot means
        if strcmp(opt.meanSymbol, '-')
            % Plot condition central tendency mean as a flat horizontal line width bar
            bx = barCenters(i1,i2)+[-0.5 0.5]*opt.barWidth;
            by = ones(1,2)*nanmean(ds{i1, i2});
            hMean = plot(bx, by, '-', 'Color', squeeze(opt.meanColors(i1,i2,:)), 'LineWidth',opt.meanLineThick);
        else
            % Plot condition central tendency mean as an independent point symbol
            hMean = plot(barCenters(i1,i2), nanmean(ds{i1, i2}), opt.meanSymbol, 'Color', squeeze(opt.meanColors(i1,i2,:)), 'MarkerSize',opt.meanDotSize);
        end

        if opt.legendToIndivs
            %set legend handles based on individual subject dots 
            handles(i1,i2) = hdot;
        else
            handles(i1,i2) = hMean;
        end
        
        % Plot error bars
        if opt.doErrorBar
            if opt.symmetricErrorBar
                % Standard Error Bars (+/- SEM)
                if strcmp(opt.errorBarType, 'line')
                    plot([1 1]*barCenters(i1,i2), nanmean(ds{i1,i2})+[-1 1]*standardError(ds{i1,i2}), '-', 'Color', squeeze(opt.errorBarColors(i1,i2,:)), 'LineWidth', opt.errorBarWidth);
                elseif strcmp(opt.errorBarType, 'box')
                    rectWid = opt.barWidth;
                    rectHei = 2*standardError(ds{i1,i2});
                    rectX = barCenters(i1,i2)-rectWid/2; % Lower left coordinate start position
                    rectY = nanmean(ds{i1,i2}-rectHei/2);
                    rectangle('Position', [rectX rectY rectWid rectHei],'EdgeColor', [0 0 0]);
                end
            else
                % Asymmetric Bootstrapped Confidence Interval Bars
                theseDats = ds{i1,i2};
                theseDats = theseDats(~isnan(theseDats));
                CI = boyntonBootstrap(@mean,theseDats,1000,opt.errorBarCI,true);
                
                if strcmp(opt.errorBarType, 'line')
                    plot([1 1]*barCenters(i1,i2), CI,'-','Color', squeeze(opt.errorBarColors(i1,i2,:)), 'LineWidth', opt.errorBarWidth);
                    % Add horizontal capping "serifs" to bounds
                    if opt.errorBarSerifs
                        for td = 1:2
                            plot(barCenters(i1, i2)+[-1 1]*0.02*diff(xlims), CI([td td]),'-','Color', squeeze(opt.errorBarColors(i1,i2,:)), 'LineWidth', opt.errorBarWidth);
                        end
                    end
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

%% Clean Ticks, Axis Titles, and Legends
l1Centers = mean(barCenters,2);
set(gca,'XTick',l1Centers,'YTick',opt.yticks);
if isfield(opt,'xTickLabs')
    set(gca,'XTickLabel',opt.xTickLabs);
end

if ~opt.doYTickLab
    set(gca,'YTickLabel',{});
elseif isfield(opt,'yTickLabs')
    set(gca,'YTickLabel',opt.yTickLabs);
end
if isfield(opt,'xLab')
    xlabel(opt.xLab);
end
if isfield(opt,'yLab') && opt.doYLab
    ylabel(opt.yLab);
end

xlim(xlims);
ylim(opt.ylims);

% Generate legend tracking if targets exist
if isfield(opt,'legendLabs') && opt.doLegend
    [legH] = legend(handles(opt.lev1ForLegend,:),opt.legendLabs,'Location',opt.legendLoc,'AutoUpdate','off');
    if ~isempty(opt.legendTitle)
        legH.Title.String = opt.legendTitle;
    end
end

end
%TEST pairedSampleDensityPlot
% Makes 4 test panels and runs a few numerical checks on the version of
% pairedSampleDensityPlot edited in Nov 2026 (CI line position, maxXs
% output, and x-limits when the midline is not at 0).
clear; close all;

rng(1); %reproducible fake data

%% Generate sample data: two conditions, N subjects each
N = 30;
data = {randn(N,1)*1.0 + 0.5, randn(N,1)*1.3 + 1.5};

%% Base options (pairedSampleDensityPlot has no defaults, so set them all)
opt.fixKernelWidth    = false; %let ksdensity pick widths, then scale them
opt.kernelWidthFactor = 1;
opt.fixedKernelWidth  = 0.3;   %only used if fixKernelWidth = true
opt.midlineX          = 0;

c1 = hsv2rgb([0.3 0.6 0.8]);
c2 = hsv2rgb([0.7 0.6 0.8]);
opt.fillColors    = [c1; c2];
opt.edgeColors    = [c1; c2]*0.7;
opt.meanColors    = [c1; c2]*0.4;
opt.fillLineWidth = 1;

opt.plotMean      = true;
opt.meanLineWidth = 2;
opt.addCI         = true;

opt.labelXVals = true;
opt.doXLabel   = true;
opt.doLegend   = true;
opt.legendLabs = {'Cond A','Cond B'};
opt.legendLoc  = 'NorthWest';

figure('Position',[100 100 1000 800]);

%% Panel 1: default case, midline at 0
subplot(2,2,1); hold on;
[kw1, maxXs1] = pairedSampleDensityPlot(data, opt);
title(sprintf('midlineX = 0, kernel width = %.3f', kw1));
ylabel('Value');

%% Panel 2: midline shifted to 3 (tests the CI-position and maxXs fixes)
% The CI lines should sit halfway between the midline (dashed) and the
% end of each mean line. Note: x-tick labels show abs(x), not distance 
% from the midline, so they are turned off here.
opt2 = opt;
opt2.midlineX   = 3;
opt2.labelXVals = false;
subplot(2,2,2); hold on;
[kw2, maxXs2] = pairedSampleDensityPlot(data, opt2);
plot([1 1]*opt2.midlineX, ylim, 'k--');
title('midlineX = 3 (CIs should be halfway to mean)');

%% Panel 3: outline only (NaN fill colors) and a fixed kernel width
opt3 = opt;
opt3.fillColors     = NaN(2,3);
opt3.fixKernelWidth = true;
opt3.fixedKernelWidth = 0.3;
opt3.fillLineWidth  = 2;
subplot(2,2,3); hold on;
kw3 = pairedSampleDensityPlot(data, opt3);
title(sprintf('Outline only, fixed kernel width = %.2f', kw3));
ylabel('Value');

%% Panel 4: all-negative data, no mean, no CI, no legend
dataNeg = {data{1}-10, data{2}-10};
opt4 = opt;
opt4.plotMean = false;
opt4.addCI    = false;
opt4.doLegend = false;
subplot(2,2,4); hold on;
pairedSampleDensityPlot(dataNeg, opt4);
title('All-negative data, no mean/CI/legend');

%% Numerical checks
fprintf(1,'\n--- pairedSampleDensityPlot checks ---\n');

% 1. maxXs should be the peak densities, and should not depend on midlineX
peakDens = zeros(1,2);
for jj = 1:2
    dens = ksdensity(data{jj},'kernel','normal','Bandwidth',kw1);
    peakDens(jj) = max(dens);
end
fprintf(1,'maxXs (midline 0):  [%.4f %.4f]\n', maxXs1);
fprintf(1,'maxXs (midline 3):  [%.4f %.4f]\n', maxXs2);
fprintf(1,'peak densities:     [%.4f %.4f]\n', peakDens);
check1 = max(abs(maxXs1-peakDens))<1e-10 && max(abs(maxXs2-peakDens))<1e-10;
fprintf(1,'Check 1 (maxXs = peak density, independent of midline): %s\n', passFail(check1));

% 2. x-limits in panel 2 should be centered on the midline
subplot(2,2,2);
xl = xlim;
check2 = abs(mean(xl) - opt2.midlineX) < 1e-10;
fprintf(1,'Check 2 (x-limits centered on midlineX=3, got [%.3f %.3f]): %s\n', xl, passFail(check2));

% 3. CI lines in panel 2 should be halfway between midline and mean-line ends
ax = gca;
lns = findobj(ax,'Type','line');
ciXs = [];
for li = 1:numel(lns)
    xd = lns(li).XData;
    if numel(xd)==2 && xd(1)==xd(2) && xd(1)~=opt2.midlineX %vertical, not the midline
        ciXs(end+1) = xd(1); %#ok<SAGROW>
    end
end
meanEnds = [];
for li = 1:numel(lns)
    xd = lns(li).XData; yd = lns(li).YData;
    if numel(xd)==2 && yd(1)==yd(2) && xd(1)==opt2.midlineX %horizontal mean lines start at midline
        meanEnds(end+1) = xd(2); %#ok<SAGROW>
    end
end
expectedCIx = opt2.midlineX + 0.5*(meanEnds - opt2.midlineX);
check3 = numel(ciXs)==2 && numel(meanEnds)==2 && ...
    all(abs(sort(ciXs) - sort(expectedCIx)) < 1e-10);
fprintf(1,'Check 3 (CI lines halfway between midline and mean): %s\n', passFail(check3));
fprintf(1,'\n');

function s = passFail(tf)
if tf, s = 'PASS'; else, s = 'FAIL'; end
end

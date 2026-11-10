%TEST makeCorrelationMatrixFigure
% Makes 4 test figures and runs numerical checks on the version of
% makeCorrelationMatrixFigure edited in Nov 2026:
%  - y-axis ticks are right for non-square matrices (YTick = 1:nY)
%  - NaN cells are white, and real correlations near minR are NOT white
%
% Each test is in its own figure window, because the function sets the
% colormap with colormap(myColrs), which applies to the whole figure. 
clear; close all;

rng(2); %reproducible fake data

%% Test 1: square matrix, NaN diagonal, one cell exactly at minR (two-hue colormap)
nV = 5;
D = randn(40, nV);
D(:,2) = D(:,1) + 0.5*randn(40,1);   %make some correlations sizeable
D(:,4) = -D(:,3) + 0.7*randn(40,1);
[rhos1, pvals1] = corr(D);
rhos1(logical(eye(nV))) = NaN;       %blank out the diagonal
maxR = 0.8; minR = -0.8;
rhos1(5,1) = minR; rhos1(1,5) = minR; %a real correlation right at the bottom of the color scale
labs1 = {'A','B','C','D','E'};

figure('Name','Test 1: square, NaN diagonal');
makeCorrelationMatrixFigure(rhos1, pvals1, labs1, labs1, true, maxR, minR);
title('Square, NaN diagonal; cells (5,1),(1,5) = minR', 'FontWeight','normal');
res(1) = checkCorrFig('Test 1', rhos1, labs1, labs1);

%% Test 2: non-square matrix (3 rows x 6 columns), no NaNs
A = randn(40,3); B = [A(:,1)+randn(40,1), randn(40,5)];
[rhos2, pvals2] = corr(A, B);       %3 x 6
xLabs2 = {'x1','x2','x3','x4','x5','x6'};
yLabs2 = {'y1','y2','y3'};

figure('Name','Test 2: non-square');
makeCorrelationMatrixFigure(rhos2, pvals2, xLabs2, yLabs2, true);
title('Non-square (3 x 6): all 3 y labels should appear', 'FontWeight','normal');
res(2) = checkCorrFig('Test 2', rhos2, xLabs2, yLabs2);

%% Test 3: all-positive, single-hue colormap, upper triangle NaN, one cell at minR = 0
rhos3 = abs(rhos1);
rhos3(triu(true(nV))) = NaN;
rhos3(4,2) = 0;                      %exactly at minR
figure('Name','Test 3: single hue');
makeCorrelationMatrixFigure(rhos3, pvals1, labs1, labs1, true, 1, 0);
title('Single hue, minR = 0; cell (4,2) = 0 should NOT be white', 'FontWeight','normal');
res(3) = checkCorrFig('Test 3', rhos3, labs1, labs1);

%% Test 4: same as test 3 but grayscale
figure('Name','Test 4: achromatic');
makeCorrelationMatrixFigure(rhos3, pvals1, labs1, labs1, true, 1, 0, true);
title('Grayscale; NaN cells white, others gray', 'FontWeight','normal');
res(4) = checkCorrFig('Test 4', rhos3, labs1, labs1);

%% Test 5: data extend below the minR you ask for (minR = -0.5, some r < -0.5), with NaN diagonal
% Real correlations below minR should be drawn in the lowest real color, not white.
figure('Name','Test 5: data below minR');
makeCorrelationMatrixFigure(rhos1, pvals1, labs1, labs1, true, 0.5, -0.5);
title('minR = -0.5 but some r < -0.5: those should NOT be white', 'FontWeight','normal');
res(5) = checkCorrFig('Test 5', rhos1, labs1, labs1);

%% Summary
fprintf(1,'\n=== %i of %i tests passed all checks ===\n\n', sum(res), numel(res));


%% Helper: check ticks and colors in the current axes
function allPass = checkCorrFig(name, rhos, xLabs, yLabs)
ax = gca;
nX = numel(xLabs); nY = numel(yLabs);
fprintf(1,'\n--- %s ---\n', name);

% 1. tick positions
c1 = isequal(ax.XTick, 1:nX) && isequal(ax.YTick, 1:nY);
fprintf(1,'Check 1 (XTick = 1:%i, YTick = 1:%i): %s\n', nX, nY, passFail(c1));

% 2 & 3. which colormap row each cell is drawn with
img  = findobj(ax,'Type','image');
C    = img.CData;
cl   = ax.CLim;
cmap = colormap(ax);
m    = size(cmap,1);
idx  = fix((C - cl(1))/(cl(2)-cl(1))*m) + 1;  %MATLAB's scaled color mapping
idx  = min(max(idx,1),m);
isWhite = all(reshape(cmap(idx(:),:),[],3)==1, 2);
isWhite = reshape(isWhite, size(C));

nanCells = isnan(rhos);
c2 = all(isWhite(nanCells));
fprintf(1,'Check 2 (all %i NaN cells drawn white): %s\n', nnz(nanCells), passFail(c2));

c3 = ~any(isWhite(~nanCells));
fprintf(1,'Check 3 (no real correlation drawn white): %s\n', passFail(c3));

% 4. one text label per non-NaN cell
txts = findobj(ax,'Type','text');
txts = txts(~strcmp({txts.String}, '')); %ignore any empty text objects
c4 = numel(txts) == nnz(~nanCells);
fprintf(1,'Check 4 (%i value labels for %i non-NaN cells): %s\n', numel(txts), nnz(~nanCells), passFail(c4));

allPass = c1 && c2 && c3 && c4;
end

function s = passFail(tf)
if tf, s = 'PASS'; else, s = 'FAIL'; end
end

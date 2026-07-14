function makeCorrelationMatrixFigure(rhos, pvals, xLabs, yLabs, addColorBar, maxR, minR, achromatic)

if ~exist('maxR', 'var')
    maxR = max(abs(rhos(:)));
end
%if minR not provided, assume its negative maxR

if ~exist('minR', 'var')
    minR = -maxR;
end

if ~exist('achromatic','var')
    achromatic = false;
end

doTwoHues = minR<0 && maxR>0;

 
%% setup colormap
if doTwoHues
    twoHues = [179 348]/360;
    Ns = [100 101];
    hues = [ones(Ns(1), 1)*twoHues(1); ones(Ns(2),1)*twoHues(2)];
    sats = [linspace(1,1/Ns(1),Ns(1)) linspace(0, 1, Ns(2))]';
    vals = ones(size(sats))*0.7;

    myColrs = hsv2rgb([hues sats vals]);
    achromatic = false;

else
    Ns = 200;

    if achromatic %grayscale! 
        hues = ones(Ns, 1);
        sats = zeros(Ns, 1);
        vals = linspace(0.8,0.05,Ns)';

    else
        hues = 348/360*ones(Ns, 1);
        sats = linspace(0, 1, Ns)';
        vals = ones(size(sats))*0.7;
    end
    myColrs = hsv2rgb([hues sats vals]);
end


% % add a special color for NaN correlations 
nanColr = [1 1 1];
nanRho = -maxR-0.1;


emptyCells = isnan(rhos);
if any(emptyCells(:))
    rhos(emptyCells) = nanRho;
    myColrs = [nanColr; myColrs];
end

if ~achromatic
    lineColr = [0 0 0];
else
    lineColr = [1 1 1];
end




%% make figure -- oddly typing "hold on" before the imagesc command flips the dimensions! 

imagesc(rhos,[minR maxR]);
colormap(myColrs);
if addColorBar, colorbar; end

hold on;

nX = length(xLabs);
nY = length(yLabs);

%add digits
for xi = 1:nX
    %plot vert dividing line
    plot([xi xi]+0.5, [0.5 nY+0.5], '-', 'Color', lineColr);
    for yi = 1:nY
        if ~emptyCells(yi, xi)
            if pvals(yi,xi)<0.01
                sigTxt = '**';
            elseif pvals(yi,xi)<0.05
                sigTxt = '*';
            else
                sigTxt = '';
            end
            text(xi,yi,sprintf('%.2f%s',rhos(yi,xi),sigTxt),'Color',[1 1 1],'HorizontalAlignment','Center','VerticalAlignment','Middle');
        end
        %plot horiz dividing line
        plot([0.5 nX+0.5], [yi yi]+0.5, '-', 'Color', lineColr);
    end
end

%starting lines at very edge
plot([0.5 0.5], [0.5 nY+0.5], '-', 'Color', lineColr);
plot([0.5 nX+0.5], [0.5 0.5], '-', 'Color', lineColr);


xlim([0.5 nX+0.5]);
ylim([0.5 nY+0.5]);

set(gca, 'XTick', 1:nX); % center x-axis ticks on bins
set(gca, 'YTick', 1:nX); % center y-axis ticks on bins
set(gca, 'XTickLabel', xLabs,'XAxisLocation','Top'); % set x-axis labels
%xtickangle(45);
set(gca, 'YTickLabel', yLabs); % set y-axis labels


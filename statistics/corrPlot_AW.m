%% function corrPlot_AW(D, vars, doPartial)
% Makes a plot of correlation matrix, like matlab's corrplot function. 
% Inputs: 
% - D: matrix of data like you might but into the "corr" function, with
%  observations in rows and variables in columns. This function computes and
%  plots pairwise correlations between all pairs of columns. 
% 
% - vars: cell array containing the names of each column
% 
% - doPartial: whether to use partialcorr (if true) or corr function (if
% false). partialcorr returns the sample linear partial correlation coefficients between 
% pairs of variables in D, controlling for the remaining variables in D.  
% 
% - nanRowOption: the "rows" variable for corr,  If "all" then, corr tries
% to use all rows. If "complete", use only rows of the input with no
% missing values; if "pairwise", compute rho(i,j) using rows with no missing values in column i or j.

function [rhos, pvals] = corrPlot_AW(D, vars, doPartial, nanRowOption)

nV = size(D,2);


if ~exist('nanRowOption', 'var')
    nanRowOption = 'complete';
end


if ~doPartial
    [rhos, pvals] = corr(D, 'Rows', nanRowOption);
else
    [rhos, pvals] = partialcorr(D, 'Rows', nanRowOption);
end
figure; 

for ri = 1:nV
    for ci=1:nV
        subi=sub2ind([nV nV],ci,ri);
        subplot(nV, nV,subi); hold on;
        
        if ri==ci
            %histogram on the diagonal 
            histogram(D(:,ri),20);

        elseif ri>ci
            xs = D(:,ci);
            ys = D(:,ri);

            plot(xs,ys,'.');

            xrng = get(gca,'XLim');
            yrng = get(gca,'YLim');
            
            xrng = xrng+[-1 1]*0.1*diff(xrng);
            yrng = yrng+[-1 2]*0.1*diff(yrng);
            set(gca,'XLim',xrng,'YLim',yrng);

            textx = xrng(1)+0.05*diff(xrng);
            texty = yrng(1)+0.88*diff(yrng);

            if pvals(ri,ci)<0.001
                star = '***';
            elseif pvals(ri,ci)<0.01
                star = '**';
            elseif pvals(ri,ci)<0.05
                star = '*';
            else
                star = '';
            end

            %text(textx, texty, sprintf('r=%.3f%s', rhos(ri,ci), star));
            title(sprintf('r=%.3f%s', rhos(ri,ci), star))
        else
            axis off;
        end

        if ri==nV
            xlabel(vars{ci});
        end
        if ci==1
            ylabel(vars{ri});
        end

    end

end

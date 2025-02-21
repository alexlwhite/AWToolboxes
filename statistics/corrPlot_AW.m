function corrPlot_AW(D, vars)

nV = size(D,2);

[rhos, pvals] = corr(D);

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

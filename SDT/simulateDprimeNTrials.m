trueDs = [1 2 3]; %true dprime;
nd = length(trueDs);

trueBias =  0; %bias here is the distance of the criterion from the neutral point, d'/2

doCorrections = [false true]; %whether to correct for HR or FR at 0 or 1. 
                               %when doCorrection is false, we just exclude
                               %those simulations from the average. 

%numbers of trials to simulate                               
nts = [10:2:300];
nn = length(nts);

%number of simulations per num. trials
nsim = 10000;

figure; hold on;
nRows = 2;
nCols = 3;
colrs = {'r','g','b'};

%loop thru whether or not to correct HR and FR when they are 0 or 1
for dci = 1:2
    doCorrection  = doCorrections(dci);
    rowNum = dci;

    %loop thru true dprime levels
    for di=1:nd
        trueD = trueDs(di); %true d-prime
        trueC = trueD/2 + trueBias; %true criterion, relative to 0 (where 0 is the mean of the target-absent distribution)

        %true HR and FR, defined by dprime and criterion
        trueFR = 1-normcdf(trueC, 0, 1);
        trueHR = 1-normcdf(trueC, trueD, 1);

        %simulate estimates of hit rate and false alarm rate for each number of trials
        hrs = NaN(nn, nsim);
        frs = NaN(nn, nsim);
        for ni = 1:nn
            thisN = nts(ni)/2; %half the number of trials for each stim category
            hrs(ni, :) = binornd(thisN, trueHR, 1, nsim)/thisN;
            frs(ni, :) = binornd(thisN, trueFR, 1, nsim)/thisN;

            if doCorrection %correct for HR or FR equal to 0 or 1 
                %assume that had we run 2x as many trials, we would have got 1
                %response different from what we got
                rateCorr = 1/(2*thisN);
                hrs(ni, hrs(ni, :)==1) = 1-rateCorr;
                hrs(ni, hrs(ni, :)==0) = rateCorr;

                frs(ni, frs(ni, :)==1) = 1-rateCorr;
                frs(ni, frs(ni, :)==0) = rateCorr;
            end
        end

        %estimate dprime
        ds = norminv(hrs) - norminv(frs);

        %estimate criterion (relative to 0)
        cs = norminv(1-frs);

        %estimate bias: distance of c from the "neutral point", d'/2
        biases = abs(cs-ds/2);

        %now deal with estimation problems when either the hit rate or false
        %alarm rate is 0 or 1, making d' end up as inf
        isBad = hrs==1 | hrs==0 | frs==1 | frs==0;
        pBad = mean(isBad, 2);

        %If we didn't correct, exclude any dprimes with hit rate or false alarm rates are 0 or 1 (d ends up as inf)
        if ~doCorrection
            ds(isBad) = NaN;
            cs(isBad) = NaN;
            biases(isBad) = NaN;
        end

        %average over simulations:
        md = mean(ds, 2, 'omitnan');
        mb = mean(biases, 2, 'omitnan');

        %Plot mean d' estimate
        subi = sub2ind([nCols nRows], 1, rowNum);
        subplot(nRows, nCols, subi); hold on;
        plot([min(nts) max(nts)], [trueD trueD], ':', 'Color', colrs{di});
        plot(nts, md, '-', 'Color', colrs{di}, 'LineWidth',1.5);
        xlabel('N Trials');
        ylabel('Mean dprime');
        ylim([0 max(trueDs)+0.5]);
        if ~doCorrection
            title('d'' (uncorrected)');
        else
            title('d'' (corrected)');
        end

        %Plot mean bias estimate
        subi = sub2ind([nCols nRows], 2, rowNum);
        subplot(nRows, nCols, subi); hold on;
        plot([min(nts) max(nts)], [trueBias trueBias], ':', 'Color', colrs{di});
        plot(nts, mb, '-', 'Color', colrs{di}, 'LineWidth',1.5);
        xlabel('N Trials');
        ylabel('Mean |bias|');
        ylim([-0.1 0.35]);

        if ~doCorrection
            title('|bias| (uncorrected)');
        else
            title('|bias| (corrected)');
        end
        
        %plot proportion of simulations with unuseable estimates of HR and
        %HR (only applies of not corrected) 
        subi = sub2ind([nCols nRows], 3, rowNum);
        subplot(nRows, nCols, subi); hold on;
        plot(nts, pBad, '-', 'Color', colrs{di}, 'LineWidth',1.5);
        ylim([-0.1 1]);

        xlabel('N Trials');
        ylabel('p(inf d'')');
        title('p(trials) unuseable estimates')
    end
end
%% ok so with a SMALL number of trials dprime estimates are low, with a moderate number of trials d' estimates are high.  more so with fewer trials.
%% the problem is noisy estimation of HR and FR.
% If HR and FR are both underestimated or both overestimated, then that
% shows up mostly as an criterion shift, in your estimate. that will happen
% half the time.
% But if HR and FR are misestimated in opposite directions, that will show
% up as an *over*estimate of dprime. That will happen roughly half the
% time, and therefore on average a slight increase in dprime estimates.
% The fewer trials, the greater the variance in estimates of HR and FR, so
% the bigger this problem is.
%
%% absolute value of bias is always overestimated (but avergages to 0 always if don't take abs value)

%% There should be a way to correct for it
function [AccCorr, conditionalAgs, corrsValsByIndex, AgsValsByIndex] = computeAccuracyContingencies(resp,pres,nRespLevs)

%Compute correlations in accuracy for dual-task responses, depending on
%target presence on the two sides 

respPres = resp>(nRespLevs/2);
respCorrect = respPres == pres;

nT = size(resp,1);

condLabs = {'any','neither','side1Only','side2Only','both'};

nConds = length(condLabs);
AccCorr = NaN(1,nConds);
for targPres = 1:nConds
    switch condLabs{targPres}
        case 'any'
            subTrials = true(nT,1);
        case 'neither'
            subTrials = ~pres(:,1) & ~pres(:,2); 
        case 'side1Only'
            subTrials = pres(:,1) & ~pres(:,2); 
        case 'side2Only'
            subTrials = ~pres(:,1) & pres(:,2);
        case 'both'
            subTrials = pres(:,1) & pres(:,2);
    end
    
    rhos = corr(respCorrect(subTrials, :));
    AccCorr(targPres) = rhos(2);
end

corrsValsByIndex.targetsPresent = condLabs;

%Compute Ag depending on whether the other side was responded to correctly
%or incorrectly 
sides = [NaN 1 2];
otherSidesPres = [NaN 0 1];
otherCorrects = [0 1];

conditionalAgs = NaN(length(sides),length(otherSidesPres),length(otherCorrects));

for sideI=1:length(sides)
    side = sides(sideI);

    if ~isnan(side)
        thisSidePres = pres(:,side); 
        otherSidePresVals = pres(:,3-side);
        
        thisSideResp = resp(:,side); 
        otherSideCorr = respCorrect(:,3-side);
    else
        thisSidePres = [pres(:,1); pres(:,2)]; 
        otherSidePresVals = [pres(:,2); pres(:,1)];
        
        thisSideResp = [resp(:,1); resp(:,2)];
        otherSideCorr = [respCorrect(:,2); respCorrect(:,1)];
    end
    
        
    for otherPresI = 1:length(otherSidesPres)
        
        otherSidePres = otherSidesPres(otherPresI);
        
        if isnan(otherSidePres)
            otherPTrials = true(size(thisSideResp));
        else
            otherPTrials = otherSidePresVals==otherSidePres;
        end
        
        for otherCorrI = 1:length(otherCorrects)
            otherCorr = otherCorrects(otherCorrI);
            if isnan(otherCorr)
                otherCorrTrials = true(size(thisSideResp));
            else
                otherCorrTrials = otherSideCorr==otherCorr;
            end
            
            subTrials = otherPTrials & otherCorrTrials;
            
            [hr,fr] = computeROCRates(thisSidePres(subTrials),thisSideResp(subTrials),1:nRespLevs);
            conditionalAgs(sideI,otherPresI,otherCorrI) = computeAROC(hr,fr);
            
            
        end
    end
end
        
AgsValsByIndex.sides = sides;
AgsValsByIndex.otherSidesPres = otherSidesPres;
AgsValsByIndex.otherCorrects = otherCorrects;

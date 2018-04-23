%% function [betas, rSqrs, conds, regressors, badVox] = runGLM(data, PRT, conds, hrf, doLinPred, motionParams, motionParamNames, censorTRs_HighFD, avgVoxPSCBeforeGLM, doPlot, voiName) 
% run a basic GLM on MRI data 
% 
% Inputs: 
% - data: a TxV matrix of raw MRI data, where T is the number of TRs and V is
% the number of voxels. This gets converted to percent signal change 
% - PRT: a structure containing the stimulus protocol. This should contain two fields: 
%     (1) eventSequence: a Tx1 vector of integers that describes which conditions
%         were on at which time points.
%     (2) condLabs: a cell array of condition labels. 
% - conds: a vector of integers specifying which event types to pull out of eventSequence. 
%   Simple case: if you have 3 conditions in eventSequence, then conds = 1:3. 
% - hrf: a structure containing the hemodynamic response function. 
%   At minimum, it should contain two fields: 
%     (1) t, a vector of time points in seconds. 
%     (2) h, a vector of hemodynamic response values equal in size to t. 
% - doLinPred: a Boolean variable determining whether to add a
%   linear trend nuisance predictor to the regression.  
% - motionParams: a TxP matrix of head motion parameters to be included as nuissance regressors. 
%   T is the number of TRs and P is the number of motion paramteres to be included in
%   the GLM as regressors. Input an empty matrix here if no motion nuissance regressors needed.
% - motionParamNames: a 1XP cell array of names of motion parameters
% - censorTRs_HighFD: boolean, whether to "censor" TRs with FD over 0.9,
%   following Siegel et al, HBM 2014
% - avgVoxPSCBeforeGLM: a Boolean variable determing whether to average
%   percent signal change over voxels bere running GLM. Set to false if you
%   want beta weights for individual voxels in the output, true if you're
%   going to average over voxels anyway. "bad" voxels with all 0 data are
%   excluded from this mean. 
% - doPlot: a variable determining whether to plot the design
%   matrix, predictors, and best-fitting voxel timecourse, each in a new
%   figure window. doPlot can be true or false, and it serves as an integer > 0
%   to specify the beta weight used to pick the "best-fitting" voxel. So if
%   doPlot = 2, the function will plot the response of an example voxel that has the highest 
%   beta for regression number 2. 
% - voiName: name of the VOI file used to extract current data. Just used
%   to add a title to a plot, if doPlot. 
% 
% Outputs: 
% - betas: a VxN matrix of beta weights, where V is the number of voxels
%   and N is the number of predictors. 
% - rSqrs: a Vx1 vector of r-squared values, reflecting quality of model
%   fits. 
% - conds: the vector of condition indices from eventSequence that were
%   actually used. This is the intersection of the input "conds" and the
%   values that are actually in eventSequence. So, some conditions you
%   requested may not have been included in the GLM (and are therefore not in
%   betas) because they werent included in the scan. 
% - regressors = a 1xN cell array of character strings indicating the
%   names of each predictor. 
% - badVox: vector of indices of "bad" voxels with timecourses that are all
%   zeros

function [betas, rSqrs, conds, regressors, badVox] = runGLM(data, PRT, conds, hrf, doLinPred, motionParams, motionParamNames, censorTRs_HighFD, avgVoxPSCBeforeGLM, doPlot, voiName) 

%Make sure that we only use conds that are actually in the PRT. Otherwise
%some design matrix columns will be all 0s, which makes no sense. 
usedConds = unique(PRT.eventSequence); 
conds = intersect(conds, usedConds); 

nCond = length(conds);
numTRs = size(data,1);
numVox = size(data,2);

%if 0 is included as conds, that should be the 'rest' or 'blank' condition.
%Make sure it's at the end 
if any(conds==0)
    conds = conds(conds~=0);
    conds = [conds 0];
end

%% Convert data to percent signal change 
PSC = fMRtoPSC(data);

%note bad voxels
badVox = [];
for vi=1:numVox
    if all(data(:,vi)==0)
        badVox = [badVox vi];
    end
end
if ~isempty(badVox)
    fprintf(1,'\n(runGLM) Note: %i bad voxels with all 0s. Voxel indices:\n',length(badVox));
    fprintf(1,'\t%i',badVox);
    fprintf(1,'\n');
end
goodVox = setdiff(1:size(PSC,2),badVox);

%% average percent signal change over voxels, if requested  (and exclude bad voxels) 
if avgVoxPSCBeforeGLM
   PSC = mean(PSC(:,goodVox),2);
   numVox = 1;
end

%% Make HRF (HDR)

%HRF should be loaded in already, in structure hrf
lengthHRF = length(hrf.t);

% We can recover estimates of our four amplitudes by multiplying our design
% matrix a matrix with nCond columns, each containing a shifted, unscaled
% estimate of the hdr:
hMat = zeros(lengthHRF*nCond,nCond);
id = 1:lengthHRF;
for i=1:nCond
    hMat(id,i) = hrf.h';
    id = id+lengthHRF;
end

%% Make design matrix 
lenPRT = length(PRT.eventSequence);
if lenPRT~=numTRs
    lenDiff = numTRs - lenPRT;
    if lenDiff>0
        fprintf(1,'\n(runGLM) Warning: PRT.eventSequence is %i TRs shorter than actual data set.',lenDiff);
        fprintf(1,'\n(runGLM) Adding %i blank TRs (0s) to eventSequence\n',lenDiff);
        PRT.eventSequence = [PRT.eventSequence zeros(1,lenDiff)];
    else
        fprintf(1,'\n(runGLM) Warning: PRT.eventSequence is %i TRs LONGER than actual data set',abs(lenDiff));
        if all(PRT.eventSequence((lenPRT+lenDiff+1):lenPRT)==0)
            PRT.eventSequence = PRT.eventSequence(1:numTRs); 
            fprintf(1,'\n(runGLM) But those extra TRs are luckily all blanks, so just trimming them off of eventSequence.\n\n');
        else
            keyboard
        end
    end
end
plotDM = false;
X = makeDesignMatrix(PRT.eventSequence', conds, lengthHRF, plotDM); 

%if 0 is in "conds", then the blank/rest condition is added as a regressor,
%at the end of the list
zeroI = find(conds==0); 
if ~isempty(zeroI)
    blankCondI = find(strcmp(PRT.condLabs,'blank') | strcmp(PRT.condLabs,'rest') | strcmp(PRT.condLabs,'break'));
    if length(blankCondI)>1 || isempty(blankCondI)
        error('(runGLM) WARNING: unclear what label to give condition 0 as requested in conds\n');
    else
        conds(zeroI) = blankCondI;
    end
end
    
regressors = PRT.condLabs(conds);

%% Make predictors

% First, convolve design matrix with HRF: 
XX = X*hMat;

% Then rescale so that predictors go up to 1  
XX = XX./repmat(max(XX),size(XX,1),1);

% Add linear predictor? (not if data already high-pass temporal filtered) 
if doLinPred
    XX = [XX,linspace(0,1,numTRs)'];
    regressors = cat(2,regressors,{'linear'});
end

% Then add the DC predictor
XX = [XX,ones(numTRs,1)];
regressors = cat(2,regressors,{'DC'});

%Then add motion param regressors: 
if ~isempty(motionParams)
   %Normalize motion params to have a max of 1
   motionParamsNorm = motionParams./repmat(max(abs(motionParams),[],1),size(motionParams,1),1);
   XX = [XX, motionParamsNorm];

   regressors = cat(2,regressors,motionParamNames);
end

%Then add mask to "censor" TRs with high FD 
if censorTRs_HighFD && ~isempty(motionParams)
   FDThresh = 0.9; 
   
   fdCol = find(strcmp(motionParamNames,'FD')); 
   if ~isempty(fdCol)
      badTRs = motionParams(:,fdCol)>FDThresh;
      XX = [XX, badTRs];
      regressors = cat(2,regressors,{'FDoverThresh'});
   end
    
end
% Each column of this matrix XX contains the unscaled
% predicted response to each isolated stimulus condition (or nuisance
% predictor)
if doPlot & ~doPlot %never mind lets not do this 
    figure; hold on
    plot(1:numTRs,XX,'LineWidth',2);
    xlabel('TR');
    ylabel('PSC');
    title('Predictors');
    legend(regressors);
    keyboard
end

%% Do regression 

betas = NaN(numVox,size(XX,2));
rSqrs = NaN(numVox,1);
%rSqrs2 = NaN(numVox,1);

% pinv(XX)*y will give us the scale values for the best linear combination of
% the columns of XX to fit the data, y.

for vi=1:numVox
    vdat = squeeze(PSC(:,vi));
    betas(vi,:) = pinv(XX)*vdat;
    
    %compute R2, reflective of model fit quality 
    predTC = XX*betas(vi,:)';
    rSqrs(vi) = corr(vdat,predTC)^2;
    
    if mod(vi,round(numVox/10))==0
        pDone = round(100*vi/numVox); 
        fprintf('(runGLM) Computed regression for %i%% of voxels\n',pDone);
    end
    
    %just to check, compute it another way by hand
%     SSTot = sum((vdat-mean(vdat)).^2);
%     SSRes =  sum((vdat - predTC).^2);
%     
%     rSqrs2(vi) = 1-SSRes/SSTot;
end

%do the two measures of r^2 agree? 
%meanR2Diff = mean(abs(rSqrs-rSqrs2)); %yes, diffs are very small

%% Plot 
if doPlot ~= 0 
    figure; hold on;

    %A negative value of doPlot causes a plot of the voxel with the highest
    %beta weight for the regressor indexed by the absolute value of doPlot
    if doPlot<0
        b1s = betas(:,abs(doPlot)); %pick voxel based on beta weight for the regressor indexed by doPlot
        maxVox = find(abs(b1s)==max(abs(b1s)));
    
        voxTC = PSC(:,maxVox);
        predTC = XX*betas(maxVox,:)';
    
        plot(1:numTRs,predTC,'b.-');
        plot(1:numTRs,voxTC,'g.-');
        legend({'predicted','data'});
    %if doPlot == 1, then plot the mean response across voxels and the mean
    %prediction 
    elseif doPlot==1 %plot mean 
        meanTC = nanmean(PSC,2); 
        predTC = XX*nanmean(betas,1)';
       
        plot(1:numTRs,predTC,'b.-');
        plot(1:numTRs,meanTC,'g.-');
        legend({'predicted mean','mean data'});        
    end
    
    xlabel('TRs');
    ylabel('PSC');
    voiName(voiName=='_')= ' ';
    title(voiName);
        
end

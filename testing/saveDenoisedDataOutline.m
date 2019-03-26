%% code bits for doing GLMdenoise and saving cleaned VTCs

%% parameters for glmDenoise

TR          = 2;       % length of TR in seconds
stimDur     = 4;  % length of stimulus block in seconds

%whether to have glmDenoise try to optimize shape of HRF (each time-point
%is a free parameter)
optimizeHRF = false; 

if optimizeHRF
    hrfmodel = 'optimize';
    hrfknobs = [];
else
    hrfmodel = 'assume';
    hrfknobs = getcanonicalhrf(stimDur,TR)';
end

opt.numboots       = 200; %number of bootstrapping
opt.wantparametric = 1; %whether to get back the final design matrix and beta weights (not in %signal change?) in results.parametric
opt.numpcstotry    = 15; %max number of principal components for noise regressors
opt.denoisespec    = '10101'; %remove polynomial drift as well as noise


    %         <denoisespec> (optional) is a binary string or cell vector of binary strings
    %       indicating the components of the data to return in <denoiseddata>.  The
    %       format of each string should be 'ABCDE' where A indicates whether to include
    %       the signal (estimated hemodynamic responses evoked by the experiment), B
    %       indicates whether to include the polynomial drift, C indicates whether
    %       to include any extra regressors provided by the user, D indicates
    %       whether to include the noise regressors, and E indicates whether
    %       to include the residuals of the model.  If multiple strings are provided,
    %       then separate copies of the data will be returned in the rows of
    %       <denoiseddata>.  Default: '11101' which indicates that all components of
    %       the data will be returned except for the component corresponding to the
    %       estimate of the contribution of the noise regressors.



%set the directory where glmDenoise will store data 
denoiseFigDir = pwd; 
%% load the data for all scans 

nScans = 4; %or something

%pretend you have a list of VTC filenames
vtcFiles = cell(1,nScans);

%pre-initialize cell array to save the raw data in a cell array called "dada"
data = cell(1,nScans);

%pre-initialize cell array to save the design matrices 
design = cell(1,nScans); 

for i = 1:nScans
    %load vtc data 
    vtc = BVQXfile(vtcFiles{i});
    
    %VTCData is organized time x X x Y x Z
    %glmDenoise wants the dimensions to be: X x Y x Z x time
    data{i} = shiftdim(vtc.VTCData,1); %move the 1st dimension to the last 
    numTRs = vtc.NrOfVolumes;
    
    %save parameters of this VTC use when pulling out results (only need
    %this if you're trying to get glmDenoise beta weights from voxels
    %corresponding to a particular VOI). 
    if i==1
        vtcp.Resolution = vtc.Resolution;
        vtcp.VTCDataSize = size(vtc.VTCData);
        vtcp.XStart = vtc.XStart;
        vtcp.YStart = vtc.YStart;
        vtcp.ZStart = vtc.ZStart;
    end

    %clear the vtc from memory
    vtc.clearObject;

    %make the design matrix 
    
    %pretend you have also loaded "eventTRs", a list of TRs when each event of each
    %condition started 
    
    %X is the design matrix
    X = zeros(numTRs, nConds);
    for condi = 1:nConds
        X(eventTRs{condi},condi)  = 1;
    end
    %glmDenoise prefers the design matrix in sparse form
    design{i} = sparse(X);
    
end

%% run glmDenoise

[results, denoiseddata] = GLMdenoisedata(design,data,stimDur,TR,hrfmodel,hrfknobs,opt,denoiseFigDir);

%% create new VTCs that contain the denoised data
for i = 1:nScans
    %load in the original VTC file 
    vtcNew = BVQXfile(vtcFiles{i});
    
    %in the denoised data, each voxel time-course has mean 0. 
    %so now let's add the original mean back in: 
    oldData = vtcNew.VTCData;
    dMeans = mean(oldData,1);
    dMeans = repmat(dMeans, [size(oldData,1) 1 1 1]);
    
    %create new data with means added back, and this also requires shifting
    %back to the dimensionality required in teh VTC: 
    vtcNew.VTCData = shiftdim(denoiseddata{i},3) + dMeans; 
    
    %save it 
    newVTCName = 'someName.vtc';
    vtcNew.saveAs(newVTCName);
    vtcNew.clearObject;
end

%% Also useful: make MDM files for running "multi-study GLMs" in BV with less clicking 
%pretend you already have PRTs for each scan, in BV format 
%this would be the list of PRT filenames: 
PRTs = cell(1,nScans);
mdmFiles = cell(nScans,2);
for i = 1:nScans
    prt = BVQXfile(PRTs{i});
    
    %pull out condition names
    condNames = cell(1,prt.NrOfConditions);
    for prti=1:prt.NrOfConditions
        condNames{prti}=prt.Cond(prti).ConditionName{1};
    end
    %find which condition was the rest or blank to not include in GLM
    restCondI = find(strcmp(condNames,'blank'));
    
    %make design matrix (SDM)
    params.nvol = numTRs; %however many TRs were in the scan 
    params.prtr = 2000; %TR duration in ms
    params.rcond = restCondI; %which condition to remove: the rest/blank
        
    sdm = prt.CreateSDM(params);

    sdmName = 'something.sdm'; 
    sdm.saveAs(sdmName);

    mdmFiles(i,:) = {VTCs{i}, sdmName};

    prt.ClearObject;
    sdm.ClearObject;
end
    
%save the MDM
mdm = BVQXfile('new:mdm'); % create the basic mdm

mdm.NrOfStudies = size(mdmFiles,1);
mdm.XTC_RTC = mdmFiles;
mdm.PSCTransformation = 0;
mdm.zTransformation = 1;
mdm.SeparatePredictors = 0;
%mdm.SeparatePredictors flag:
%"codes whether predictors of equal name are either concatenated across all runs of all subjects (0),
% only across runs of the same subject, but separate across subjects (2),
% or fully separated across runs and subjects (1) 
mdm.RFX_GLM = 0; % random effects

mdmName = 'something.mdm';
mdm.saveAs(mdmName)

%also try this: 
%mdm.ComputeGLM;

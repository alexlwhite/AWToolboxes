%load a VOI with a given file name 
voi = BVQXfile(voiFileName);

%you'll also need a structure "vtcParams" that has the resoltion and size
%of the data matrix in the original VTC. So let's say you had a VTC, you
%could create the vtcParams structure as follows:
vtcParams.Resolution = vtc.Resolution;
vtcParams.VTCDataSize = size(vtc.VTCData);
vtcParams.XStart = vtc.XStart;
vtcParams.YStart = vtc.YStart;
vtcParams.ZStart = vtc.ZStart;

%get out the linear indices of voxels for that ROI 
linearFuncIs = getVTCFunctionalVoxelsFromVOI(vtcParams, voi);
linearFuncIs = linearFuncIs{1};

%pull out beta weights for those voxels from the results structure returned
%by glmDenoise
%for a particular regressor, in this case the 2nd one: 
regressorNum = 2; 

%first pull out all beta weights for that regressor
thisRegBetas = squeeze(res.modelmd{2}(:,:,:,regi));

%Then average over verage over voxels with the indices corresponding to the ROI 
betas(reti,stimROI,regi,di) = mean(thisRegBetas(linearFuncIs));


%% function PSC = fMRtoPSC(data)
% convert raw fMR data to percent signal change. 
% 
% Inputs: 
% - data: a TvV matrix, where T is the number of time points (TRs) and V is the number
% of voxels. Time in rows, voxels in columns. 
% 
% Ouputs: 
% - PSC: a TxV matrix of the data transformed voxelwise into percent signal
% change: voxelPSC = (voxelData - voxelMean)/voxelMean

function PSC = fMRtoPSC(data)

numTRs = size(data,1);

voxMeans = repmat(mean(data,1),[numTRs 1]); 
PSC = 100*(data - voxMeans)./voxMeans;
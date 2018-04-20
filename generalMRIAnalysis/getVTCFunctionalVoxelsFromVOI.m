%% [uniqueLinearIs, uniqueFuncCoords] = getVTCFunctionalVoxelsFromVOI(vtc, voi, vmrBBox)
% The purpose of this function is to return indices and coordinates of
% functional voxels corresponding to a VOI, so that data can be pulled out
% of VTC.VTCData. This is meant for use with BrainVoyager vtc and voi
% files, in combination with the Matlab BVQX tools. 
% 
% Inputs: 
% - vtc: a BVQX vtc object 
% - voi: a BVQX voi object 
% - vmrBBox: the structure returned by vmr.BoundingBox, for the anatomical
%   VMR to which these data are aligned. If this is not input, it is by
%   defualt set to vmrBBox.DimXYZ = [256, 256, 256]; 
% 
% Outputs: 
% - uniqueLinearIs: a 1xC cell array, for each of C VOIs within the input voi. 
%   For each ci, uniqueLinearIs{ci} is a Vx1 vector of linear indices of voxels 
%   within vtc.VTCData that fall in the VOI, in functional resolution. V is 
%   the number of *unique* functional voxels within the VOI. 
%    This can be used to pull out the data like so: 
%      data = double(vtc.VTCData(:,uniqueLinearIs{ci}))
% - uniqueFuncCoords: a 1xC cell array, for each of C VOIs within the input voi. 
%   For each ci, uniqueFuncCoords{ci} is a Vx3 matrix of functional coordinates 
%   of voxels within VTCData that fall in the VOI. 
% 
% by Alex White, 2018, at the University of Washington
% 

function [uniqueLinearIs, uniqueFuncCoords] = getVTCFunctionalVoxelsFromVOI(vtc, voi, vmrBBox)

if nargin < 3 || ~isstruct(vmrBBox) || numel(vmrBBox) ~= 1 || ~isfield(vmrBBox, 'DimXYZ')
    vmrBBox = struct('DimXYZ', [256, 256, 256]);
end

nVOIs = voi.NrOfVOIs;

uniqueLinearIs = cell(1,nVOIs);
uniqueFuncCoords = cell(1,nVOIs);
uvec     = cell(1,nVOIs);
uvecr    = cell(1,nVOIs);

for vc=1:nVOIs
    %Pull out 3D coordinates for VOI in anatomical space
    anatCoords = voi.BVCoords(vc, vmrBBox);
    %This is the same as voi.VOI(1).Voxels;
    
    % get VTC info
    vres = vtc.Resolution;
    vtcsz = size(vtc.VTCData);    
    voff = [vtc.XStart, vtc.YStart, vtc.ZStart];
    
    % Translate from anatomical to functional space
    funcCoords = round(1 + (anatCoords - repmat(voff, [size(anatCoords, 1), 1])) / vres);
    
    
    % remove bad entries that are out of range: 
    be = (funcCoords(:, 1) < 1 | funcCoords(:, 1) > vtcsz(2) | ...
        funcCoords(:, 2) < 1 | funcCoords(:, 2) > vtcsz(3) | ...
        funcCoords(:, 3) < 1 | funcCoords(:, 3) > vtcsz(4));
    
    funcCoords(be, :) = [];
    
    % functCoords is a nVox x 3 matrix, with functional coordinates xyx,
    % and nVox is the number of *anatomical* voxels, so there are lots of
    % repeats in this list of functional voxels
    
    % Go from 3D coordinates to linear indices
    voxs = sub2ind(vtcsz(2:4), funcCoords(:, 1), funcCoords(:, 2), funcCoords(:, 3));
    %so now voxs is a nVox x 1, where nVox is still the number of anatomical
    %voxels, but it has lots repeats, and each value is the linear index of
    %*functional* voxels to get pull from VTCData
    
    %Take unique voxel indices
    [uniqueLinearIs{vc}, uvec{vc}, uvecr{vc}] = unique(voxs);
    %uvec and uvecr are equivalent to what vtc.VOITimeCourse returns
    
    %Take unique coordinates
    uniqueFuncCoords{vc} = unique(funcCoords, 'rows');
     
    
    %% check to make sure all methods return the same data 
%     uniqueLinearIs2 = sub2ind(vtcsz(2:4), uniqueFuncCoords{vc}(:, 1), uniqueFuncCoords{vc}(:, 2), uniqueFuncCoords{vc}(:, 3));
%     uniqueLinearIs2 = unique(uniqueLinearIs2);
%     linearIMatch = all(uniqueLinearIs{vc} == uniqueLinearIs2)
%    
%     data1 = double(vtc.VTCData(:,uniqueLinearIs{vc})); 
%     
%     weightParam = inf;
%     [voitc, uvec2, uvecr2] = vtc.VOITimeCourse(voi, weightParam); % help on this: vtc.help('VOITimeCourse')
%     data2 = voitc{vc}(:,uvec2{vc});
% 
%     dataMatch = all(all(data1==data2))
%     uvecMatch = all(uvec2{vc}==uvec{vc})
%     uvecrMatch = all(uvecr2{vc}==uvecr{vc})
    
end

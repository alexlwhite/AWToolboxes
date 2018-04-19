%% function pOverlap = computeVOIOverlap(f1,f2)
% Computes the proportion overlap between two VOIs (ROI definition files) 
% pOverlap is the number of voxels in the intersection of the two VOIs,
% divided by the number of voxels in the union of the two VOIs. 
% 
% Inputs: 
% - f1: full file name of VOI number 1
% - f1: full file name of VOI number 2
% 
% Output: 
% - pOverlap: is the number of voxels in the intersection of the two VOIs,
% divided by the number of voxels in the union of the two VOIs. 

function pOverlap = computeVOIOverlap(f1,f2)

%the first 14 lines contain header info. textscan will put each item
%(separated by tab) into a new element, and each line contains 2 such
%elements 
nLines = 14;
numElementsPrint = nLines*2;

%add two elements for color 
colorIs = numElementsPrint:(numElementsPrint+2);
nElements = colorIs(end);


%element that is VOI name: 
nameI = 26;

fid1 = fopen(f1,'r'); 
fid2 = fopen(f2,'r'); 

frmt = '%s';

try
    C1 = textscan(fid1,frmt,nElements);

    C2 = textscan(fid2,frmt,nElements);
catch
    keyboard
end

%element that is naming convention...seems to have switched across BV
%versions?
conventionI1 = find(strcmp(C1{1},'<SUBJ>_<VOI>') | strcmp(C1{1},'<VOI>_<SUBJ>'));
conventionI2 = find(strcmp(C2{1},'<SUBJ>_<VOI>') | strcmp(C2{1},'<VOI>_<SUBJ>'));

headerRowMismatch =  conventionI1 ~= conventionI2;

convention1 = C1{1}{conventionI1};
convention2 = C2{1}{conventionI2};

conventionMismatch = ~strcmp(convention1, convention2);
if conventionMismatch
    fprintf(1,'\n(combineVOIs): Two VOIs have different naming conventions (<SUBJ>_<VOI> or <VOI>_<SUBJ>)\n');
end

areSame = strcmp(C2{1},C1{1});
valsToCheck = 1:(nameI-1);
if conventionMismatch && ~headerRowMismatch
    valsToCheck = setdiff(valsToCheck,conventionI1);
end

if ~all(areSame(valsToCheck))
    keyboard
    error('(combineVOIs): Two VOIs have different header info!');
end

%Use BVQX to load in the VOIs and extract voxel info
voi1 =  BVQXfile(f1);
vox1 = voi1.BVCoords(1);
voi2 =  BVQXfile(f2);
vox2 = voi2.BVCoords(1);

%voi.VOI.voxels has xyz coordinates of anatomical voxels (one row per voxel). Convert to linear
%indices: 
linIs1 = sub2ind(ones(1,3)*voi1.OriginalVMRFramingCubeDim, vox1(:,1),vox1(:,2),vox1(:,3)); 
linIs2 = sub2ind(ones(1,3)*voi1.OriginalVMRFramingCubeDim, vox2(:,1),vox2(:,2),vox2(:,3)); 

intersectIs = intersect(linIs1,linIs2);
unionIs = union(linIs1,linIs2);
pOverlap = numel(intersectIs)/numel(unionIs);

voi1.clearObject;
voi2.clearObject;
 
clear voi1 voi2

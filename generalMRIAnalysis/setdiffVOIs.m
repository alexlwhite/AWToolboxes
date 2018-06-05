%% function newFile = setdiffVOIs(f1,f2,newName,folder)
%
% Combine two BrainVoyager VOIs by taking the set difference of their
% voxels - that is, find the part of A that is not in B.
% Reads in two file ".voi" text files, and  writes a new file with the same header information, but a list of voxels
% that is the setdiff of the two input VOIs. Duplicate voxels are excluded. 
% 
% Inputs: 
% - f1: character string, the name of the .voi file for VOI 1 
% - f1: character string, the name of the .voi file for VOI 2 
% - newName: character string containing the name (just the name, without
%   directory) of the new combined .voi file. Should have extension .voi. 
% - folder: character string, directory to which the resulting .voi file
%   should be written
% 
% Outputs: 
% - newFile: the name (with directory) of the resulting .voi file 
% - numVox: the number of voxels in the ROI
% 
function [newFile, numVox]= setdiffVOIs(f1,f2,newName,folder)

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

diffIs = setdiff(linIs1,linIs2);


[i,j,k] = ind2sub(ones(1,3)*voi1.OriginalVMRFramingCubeDim, diffIs);
newVox = [i j k];
numVox = size(newVox,1);

%open new file 
newFile = fullfile(folder,newName);
fidNew = fopen(newFile,'w');

%write header
for ei=1:numElementsPrint
    if ei==nameI
        fprintf(fidNew,'%s\n',newName);
    elseif ei>=colorIs(1)
        fprintf(fidNew,'%s ',C1{1}{colorIs});
        fprintf(fidNew,'\n');
    else
        if mod(ei,2)==1
            fprintf(fidNew,'%s\t',C1{1}{ei});
        else
            fprintf(fidNew,'%s\n',C1{1}{ei});
        end
    end
end

%write voxels
fprintf(fidNew,'\nNrOfVoxels: %i\n',numVox);
for vi=1:numVox
    fprintf(fidNew,'%i %i %i\n',newVox(vi,1),newVox(vi,2),newVox(vi,3));
end

%write last line: 
fprintf(fidNew,'\nNrOfVOIVTCs: 0');
    
    
%% Clear BVQX objects
voi1.ClearObject;
voi2.ClearObject;
if exist('clearbvxqobjects')
    clearbvqxobjects({voi1,voi2});
end

fclose(fid1); fclose(fid2); fclose(fidNew);
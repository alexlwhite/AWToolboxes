%% function uniqueVoxels = excludeDuplicateVOIVoxels(voiFile,voxelsToExclude,newName,folder)
% For a given VOI (pointed to by voiFile), remove a set of voxels in list
% voxelsToExclude. Then writes a new VOI file, using the same header but new
% list of voxels with those in voxelsToExclude removed. 
% 
% Inputs: 
% - voiFile, full name with directory of a VOI txt file 
% - voxelsToExclude: list of voxels to exclude, as linear indices (not 3D
% coordinates). 
% - newName: name of the new VOI
% - folder: directory into which to write the new VOI
% 
% Outputs: 
% - uniqueIs: list of unique indices in the resulting VOI, as linear
% indices. 


% 
function uniqueIs = excludeDuplicateVOIVoxels(voiFile,voxelsToExclude,newName,folder)

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

fid1 = fopen(voiFile,'r'); 

frmt = '%s';
C1 = textscan(fid1,frmt,nElements);
  

%Use BVQX to load in the VOIs and extract voxel info
voi1 =  BVQXfile(voiFile);
vox1 = voi1.BVCoords(1); %xyz coordinates of anatomical voxels

%look for duplicates
linIs = unique(sub2ind(ones(1,3)*voi1.OriginalVMRFramingCubeDim, vox1(:,1),vox1(:,2),vox1(:,3))); 
nVox = length(linIs);

uniqueIs = unique(setdiff(linIs,voxelsToExclude)); 

nDuplicates = nVox - length(uniqueIs);
if nDuplicates>0
    fprintf(1,'\n(excludeDuplicateVOIVoxels): VOI %s had %i duplicate voxels (removed)\n',newName,nDuplicates); 
end

[i,j,k] = ind2sub(ones(1,3)*voi1.OriginalVMRFramingCubeDim, uniqueIs);
uniqVox = [i j k];
nUniqVox = size(uniqVox,1);

%Re-write the file 
fidNew = fopen(fullfile(folder,newName),'w');

%re-write header
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
fprintf(fidNew,'\nNrOfVoxels: %i\n',nUniqVox);
for vi=1:nUniqVox
    fprintf(fidNew,'%i %i %i\n',uniqVox(vi,1),uniqVox(vi,2),uniqVox(vi,3));
end

%write last line: 
fprintf(fidNew,'\nNrOfVOIVTCs: 0');
    
%% Clear BVQX objects
voi1.ClearObject;
if exist('clearbvxqobjects')
    clearbvqxobjects({voi1});
end
    
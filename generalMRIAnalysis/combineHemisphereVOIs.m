function combinedFile = combineHemisphereVOIs(f1,f2,newName,folder)

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

C1 = textscan(fid1,frmt,nElements);
C2 = textscan(fid2,frmt,nElements);

areSame = strcmp(C2{1},C1{1}); 
if ~all(areSame(1:(nameI-1)))
    error('(combineHemisphereVOIs): Two VOIs have different header info!');
end

%Use BVQX to load in the VOIs and extract voxel info
voi1 =  BVQXfile(f1);
vox1 = voi1.BVCoords(1);
voi2 =  BVQXfile(f2);
vox2 = voi2.BVCoords(1);

%voi.VOI.voxels has xyz coordinates of anatomical voxels
allVox = [vox1; vox2];

nVox = size(allVox,1);

%look for duplicates
linIs = sub2ind(ones(1,3)*voi1.OriginalVMRFramingCubeDim, allVox(:,1),allVox(:,2),allVox(:,3)); 

uniqIs = unique(linIs); 

nDuplicates = nVox - length(uniqIs);
fprintf(1,'\n(combineHemisphereVOIs): VOI %s had %i duplicate voxels (removed)\n',newName,nDuplicates); 

[i,j,k] = ind2sub(ones(1,3)*voi1.OriginalVMRFramingCubeDim, uniqIs);
uniqVox = [i j k];
nUniqVox = size(uniqVox,1);

%open new file 
combinedFile = fullfile(folder,newName);
fidNew = fopen(combinedFile,'w');

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
fprintf(fidNew,'\nNrOfVoxels: %i\n',nUniqVox);
for vi=1:nUniqVox
    fprintf(fidNew,'%i %i %i\n',uniqVox(vi,1),uniqVox(vi,2),uniqVox(vi,3));
end

%write last line: 
fprintf(fidNew,'\nNrOfVOIVTCs: 0');
    
    
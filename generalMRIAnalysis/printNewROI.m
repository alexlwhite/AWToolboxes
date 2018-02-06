%% function combinedFile = printNewROI(exampleFile,newName,folder,newVoxels)
% prints text file in VOI format for BrainVoyager 
% 
% Inputs: 
% - exampleFile: the name of another VOI text file from which header
% information can be loaded and copied 
% - newFolder: name of directory into which to put the new voi text file
% - newFile:  name of new text file to write out 
% - newVoxels: a Nx3 matrix containing anatomical coordinates of N voxels in the ROI
% 
% OutputS: 
% fidNew: handle of new text file. Will be -1 if failed to open file. 

function fidNew = printNewROI(exampleFile,newFolder,newName,newVoxels)

%the first 14 lines contain header info. textscan will put each item
%(separated by tab) into a new element, and each line contains 2 such
%elements 
nLines = 14;
numElementsPrint = nLines*2;

%add two for color 
colorIs = numElementsPrint:(numElementsPrint+2);
nElements = colorIs(end);

%element that is VOI name: 
nameI = 26;

fid1 = fopen(exampleFile,'r'); 

frmt = '%s';
C1 = textscan(fid1,frmt,nElements);


%open new file 
newFile = fullfile(newFolder,newName);
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
numVox = size(newVoxels,1);
fprintf(fidNew,'\nNrOfVoxels: %i\n',numVox);
for vi=1:numVox
    fprintf(fidNew,'%i %i %i\n',newVoxels(vi,1),newVoxels(vi,2),newVoxels(vi,3));
end

%write last line: 
fprintf(fidNew,'\nNrOfVOIVTCs: 0');
    
    
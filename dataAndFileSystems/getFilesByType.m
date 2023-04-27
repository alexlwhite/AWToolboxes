%% function [fileNames, dirs, dirIs] = getFilesByType(sdir,filetype)
%Finds all the files of a specific type in a directory, including
%subfolders. Searches recursively. 
% By Alex White, 2018 
% 
% Inputs: 
% - sdir: directory to search through 
% - filetype: character string of the file extension, like 'txt', 'mat', 'edf'
% 
% Outputs
% - fileNames: a 1xN cell array of the names of files found with the given
% type. This is the full file path. 
% - dirs: a 1xM cell array of the directories that contain files of the given type
%   immediately within them (that is, not within a subfolder)
% - dirIs: a 1xN vector of indices of directories in dirs that each file in
%   fileNames lives in. 


function [fileNames, dirs, dirIs] = getFilesByType(sdir,fileType)

fnames = sort(getFileNames(sdir));
dirs = {}; fileNames = {}; dirIs = [];
nDir = 0; nFiles = 0;

%how many characters in the 
typeLen = length(fileType);

%figure out which files in this directory are themselves directories, and
%do them last 
isFolder = NaN(1,numel(fnames)); 
for fi=1:numel(fnames) 
    thisf = fullfile(sdir,fnames{fi});
    isFolder(fi)= exist(thisf)==7;
end
fileOrder = [find(~isFolder) find(isFolder)];


for fi = fileOrder
    thisf = fullfile(sdir,fnames{fi});
    %if this is itself a subfolder
    if isFolder(fi)
        [newFs, newDirs, newIs] = getFilesByType(thisf,fileType);
        if numel(newFs)>0
            fileNames = cat(2,fileNames,newFs);
            dirs = cat(2,dirs,newDirs);
            dirIs = [dirIs newIs+nDir]; 
            nDir = nDir+numel(newDirs);
            nFiles = numel(fileNames);
        end
    else
        if length(thisf)>(typeLen+1)
            ftype = thisf((end-typeLen+1):end);
            if strcmp(ftype,fileType)
                nFiles = nFiles+1;
                fileNames{nFiles} = thisf;
                if isempty(dirs) %initalize dirs
                    dirs{1} = sdir;
                    nDir = 1;
                end
                dirIs = [dirIs 1];
            end
        end
    end
end

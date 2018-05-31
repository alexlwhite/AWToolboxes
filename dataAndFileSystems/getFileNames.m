function fileNames = getFileNames(dataDir)
% Get a list of files in a directory as a cell array
%
% fileNames = getFileNames(dataDir)
% 
% Input is a character string with name of directory to be searched
% Output is a cell array of filenames inside. File names don't contain full
% address, just name of that file. 

f = dir(dataDir);
c = 0;
fileNames = {};
for ii = 1:length(f)
    if ~strcmp(f(ii).name,'.') && ~strcmp(f(ii).name,'..') && ~strcmp(f(ii).name,'.DS_Store')
        c = c+1;
        fileNames{c} = f(ii).name;
    end
end



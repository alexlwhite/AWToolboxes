function fileNames = getSubfolders(dataDir)
% Get a list of directories inside a directory as a cell array
%
% fileNames = getSubfolders(dataDir)
% 
% Input is a character string with name of directory to be searched
% Output is a cell array of names of folders inside. Names don't contain full
% address, just name of that folder. 

f = dir(dataDir);
c = 0;
fileNames = {};
for ii = 1:length(f)
    if ~strcmp(f(ii).name,'.') && ~strcmp(f(ii).name,'..') && ~strcmp(f(ii).name,'.DS_Store') && f(ii).isdir
        c = c+1;
        fileNames{c} = f(ii).name;
    end
end



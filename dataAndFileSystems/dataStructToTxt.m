%% function dataStructToTxt(dat, filename)
% Takes a data structure and prints it all out to a text file. 
% Each field in the data structure, which should be a vector, gets 1 column 
% in the text file, and each trial gets a row. 
% 
% Inputs: 
% - dat: a structure. This function will try to print out each field as a
%  column in the txt file. This will only work if all fields in the
%  structure have the same length 
% - filename: the full name of the txt file to print to 
% - multipleRowVariable: if there is 1 data variable that may have
%   multiple values per trial, each trial will have n rows depending on how
%   many vales of this variable there are on that trial. Leave empty if that
%   doesn't apply. (Commonly used for the nSac varabiable, for data
%   structures in which there are multiple saccades and we print out stats
%   for each one). 


function dataStructToTxt(dat, filename, multipleRowVariable)

if nargin<3
    multipleRowVariable = '';
end

%open the text file
df = fopen(filename,'w');

%extract the names of data fields to export
fs = fieldnames(dat);

%first, print header
for f=1:numel(fs)
    fprintf(df,'%s\t',fs{f});
end
fprintf(df,'\n');



%how many trials: depends on what the trial counter is called
if isfield(dat,'trial')
    nTrials = length(dat.trial);
elseif isfield(dat,'t')
    nTrials = length(dat.t); 
else %if no explicit trial counter, just take the first field
    eval(sprintf('nTrials = length(dat.%s)',fs{1}));
end

for t=1:nTrials
    %determine how many repeated rows needed of this trial
    if isempty(multipleRowVariable)
        n = 1;
    else
        eval(sprintf('n = max([1 dat.%s(t)]);', multipleRowVariable));
    end
    for s=1:n
        for f=1:numel(fs) 
            eval(sprintf('docell = iscell(dat.%s);',fs{f})); 
            
            if docell
                eval(sprintf('fprintf(df,''%%3f\t'',dat.%s{t}(s));',fs{f}));
            else
                try
                    eval(sprintf('fprintf(df,''%%3f\t'',dat.%s(t));',fs{f}));
                catch me
                    fprintf(1,'\nERROR in dataStructToTxt, with variable %s\n',fs{f});
                    keyboard
                end
            end
        end
        fprintf(df,'\n');
    end
end
        
%% At start of experiment, set EDF name and open the Eyelink file

%Set file names: 
eyelinkFileName = 'XX_03_1.edf';  %something less than 8 characters
dataFileName = 'ExperimentFolder/data/XX_180103_01.edf'; %something informative, including the full path of where I want it to end up

%open file: 
i = Eyelink('OpenFile', eyelinkFileName);		% open data file on operater PC
if i~=0
    fprintf('Cannot create EDF file ''%s'' ', eyelinkFileName);
    Eyelink( 'Shutdown');
    return;
end

%% At end of experiment, transfer and close the file 

%Get edf file from the Eyelink computer to the experiment computer: 
status = Eyelink('ReceiveFile'); 
if status == 0
    fprintf(1,'\n\nFile transfer went pretty well\n\n');
elseif status < 0
    fprintf(1,'\n\nError occurred during file transfer\n\n');
else
    fprintf(1,'\n\nFile has been transferred (%i Bytes)\n\n',status)
end

%Note: there are more options for the ReceiveFile function, which I don't use. Here's the help text: 

% [status =] Eyelink('ReceiveFile',['filename'], ['dest'], ['dest_is_path'])
% 
%  If <src> is omitted, tracker will send last opened data file.
%  If <dest> is omitted, creates local file with source file name.
%  Else, creates file using <dest> as name.  If <dest_is_path> is supplied and
% non-zero
%  uses source file name but adds <dest> as directory path.
%  returns: file size if OK, 0 if file transfer was cancelled, negative =  error
% code


%Move the file from one place to another on the experiment computer
[success, message] = movefile(eyelinkFileName,dataFileName);

%Close the file and shut down 
Eyelink('closefile');



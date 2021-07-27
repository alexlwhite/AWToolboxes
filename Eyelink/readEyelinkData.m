function result = readEyelinkData(fname, deleteASC)
% result = readEyelinkData(fname, deleteASC)
% Read an EDF file or an ASC file.
% 1) If filename ends with '.edf' (case sensitive), run edf2asc - this needs to
%    be on the path!
% 2) Open resulting '.asc' file (or the specified '.asc' file) and read the
%    data.
% 3) Saccade events translated into [startTime endTime duration startX startY endX
%    endY ? ?]
% 4) Fixation events translated into [startTime endTime duration avX avY ?]
% 5) MSG events recorded: first word (message name) becomes a fieldname in the
%    structure, appended with _m (a string containing the message) or _t
%    (the time of the message). In addition, each trial has a cell array
%    messages, and corresponding vector of messageTimes.
% Timestamps are corrected to be 0 at start of each trial.
%
% written by Alex White as a modified version of Sanjay Manohar's
% readEDFASC.

if ~exist('deleteASC','var'), deleteASC=0; end

%% determine what type of file we're given (.edf or .asc)
if fname(end-3)=='.'
    fstem=fname(1:end-4);
else
    fstem=fname;
end

rehash path
if exist([fstem '.asc'],'file')
    fname=[fstem '.asc'];
    fprintf(1,'\n(%s) Using existing asc file %s\n', mfilename, fname);
elseif exist([fstem '.edf'], 'file')
    fname=[fstem '.edf'];
end

%% if input is an EDF file, run edf2asc
if strcmp(fname(end-3:end), '.edf')
    fprintf(1,'\n(%s) Running edf2asc on %s...', mfilename, fname);
    fname2=[fname(1:end-4) '.asc'];
    
    [err,tmp]=dos(['edf2asc ' fname], '-echo');
    rehash path
    if ~exist(fname2, 'file')  % error reading file
        fprintf(tmp); fprintf('\n error reading edf\n');
        result=struct();
        return;
    end
    disp(tmp);
else
    fname2=fname;
end

%% open the ASC file
[fid, mess] = fopen(fname2, 'r');
if(fid==-1), error(mess); end

%% read the file
line=0;
trial=0;
trialstart=0;
nancount=0;
maxnans=3000; %max 3 seconds blink recorded.
%chars allowed in messages
goodChars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890.,_; ';

try
    skipToNextLineStartingWith(fid, 'MSG');
    token='MSG';
    
    while ~feof(fid)
        if strcmp(token,'MSG')
            t = fscanf(fid, '%d',1); %message time
            msg = fscanf(fid, '%s',1); %message text
            val = fgets(fid); %remainder of line
            wholeMessage = [msg val];
            %remove nonalphanumeric characters:
            wholeMessage = wholeMessage(ismember(wholeMessage, goodChars));
            
            %ASSUMES THAT EACH TRIAL STARTS WITH A TRIALID MESSAGE, AND
            %THATS HOW THE DATA ARE DIVIDED UP
            if strcmp(msg, 'TRIALID')
                trial=trial+1;
                result(trial).gazeData=[];
                result(trial).fixation=[];
                result(trial).saccade=[];
                result(trial).blink=[];
                result(trial).messages = {};
                result(trial).messageTimes = [];
                trialstart=fscanf(fid,'%d',1);
                fprintf('\rTrial %d   ',trial);
            elseif(msg=='B')
                
                %not sure what this does
                tmp=find(val==':');
                if ~isempty(tmp)
                    [msg,tmp2,tmp2,tmp2]=sscanf(val(tmp+1:end),'%s',1);
                    val=val((tmp+tmp2):end);
                else
                    [result(trial).B tmp tmp tmp]= sscanf(val,'%d',1);
                    [result(trial).T tmp tmp tmp]= sscanf(val((tmp+2):end),'%d',1);
                    msg='BT';
                end
            end
            if trial>0
                fn = removeNonalphanumericChars(msg);
                val(uint8(val)<31)=[]; % remove control characters from message
                %store message "value" and time
                result(trial).([fn '_m']) = val;
                result(trial).([fn '_t']) = t-trialstart;
                
                %add just a cell array of the whole messages
                result(trial).messages = cat(1, result(trial).messages, {wholeMessage});
                result(trial).messageTimes = [result(trial).messageTimes; t-trialstart];
            end
            %start of this line is a timestamp
        elseif any(token(1)=='0123456789') && trial>0
            %pull out the gaze data
            gazeDat=fscanf(fid,'%g',3)';
            time=str2num(token)-trialstart ;
            if length(gazeDat)==3
                result(trial).gazeData = [result(trial).gazeData; time gazeDat] ;
                nancount=0;
            else
                if(nancount<maxnans)
                    result(trial).gazeData = [result(trial).gazeData; time NaN NaN NaN];
                elseif nancount==maxnans
                    fprintf('?blink?');
                end
                nancount=nancount+1;
            end
            fgets(fid); %gobble rest of line
        elseif strcmp(token,'EFIX') && trial>0
            tmp=fscanf(fid,'%s',1); t=fscanf(fid,'%g',6)';
            result(trial).fixation=[result(trial).fixation; t - [trialstart trialstart 0 0 0 0]];
            %fprintf('f');
        elseif strcmp(token,'ESACC') && trial>0
            tmp=fscanf(fid,'%s',1); t=fscanf(fid,'%g',9)';
            if length(t)==9
                result(trial).saccade=[result(trial).saccade; t - [trialstart trialstart 0 0 0 0 0 0 0]];
            else
                fprintf(['?']); disp(t);
            end;
            %fprintf('s');
        elseif strcmp(token, 'EBLINK') && trial>0
            tmp=fscanf(fid,'%s',1); t=fscanf(fid,'%g',3)';
            result(trial).blink  =[result(trial).blink; t - [trialstart trialstart 0]];
            %fprintf('B');
        else
            fgets(fid); %gobble
        end
        
        %get next "token"
        token=fscanf(fid, '%s',1);
    end
    if ~(prod(size(result(1).gazeData)))
        result(1)=[]; %remove single initial blank trial
    end
    if isfield(result(1), 'VOID_TRIAL_t')
        n=sum([result.VOID_TRIAL_t]>0);
        if ~automatic
            if(input(['Delete ' num2str(n) ' void trials? (1/0)']))
                result([result.VOID_TRIAL_t]>0)=[];
            end
        end
    end
    
catch
    e=lasterror;
    fprintf('%s\nin %s\nline %d\n',e.message, e.stack(1).file, e.stack(1).line);
end
fclose(fid);
%% delete ASC file
if deleteASC
    delete(fname2);
end

%rename some useful things
for trial=1:length(result)
    result(trial).timeStamp = result(trial).gazeData(:,1);
    result(trial).gazePosX = result(trial).gazeData(:,2);
    result(trial).gazePosY = result(trial).gazeData(:,3);
    result(trial).pupilSize = result(trial).gazeData(:,4);
    
    %extract some facts about recording
    modeMsg = result(trial).MODE_m;
    spaces = find(modeMsg == ' ');
    modeWords = {};
    for spi=1:length(spaces)
        if spi<length(spaces)
            modeWord = modeMsg((spaces(spi)+1):(spaces(spi+1)-1));
        else
            modeWord = modeMsg((spaces(spi)+1):end);
        end
        modeWords = cat(2, modeWords, {modeWord});
    end
    if ~isempty(modeWords)
        result(trial).sampleRate = modeWords{3};
        result(trial).trackedEye = modeWords{end};
    end
end
function skipToNextLineStartingWith(fid, str)
token = ''; pass=1;
while ~feof(fid) & pass
    token=fscanf(fid, '%s',1);
    if strcmp(token,str) pass=0; break;
    else fgets(fid);  %gobble
    end;
end;

function str=removeNonalphanumericChars(str)
% Remove all characters that are not alphanumeric, . or _,  and remove
% any initial digits.
% useful for making field names or variable names.
i=1;
while i<=length(str)
    if any(str(i)=='_.') i=i+1;continue; end;
    if str(i)<65 | (str(i)>90 & str(i)<97) | str(i)>122 ...
            | (i==1 & str(i)>47 & str(i)<58)
        str=[ str(1:i-1)  str(i+1:end) ];
    else i=i+1;
    end
end


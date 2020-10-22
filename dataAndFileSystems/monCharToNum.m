function monNum = monCharToNum(monChar)

months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
monNum = find(strcmp(monChar, months));

if isempty(monNum)
    fprintf(1,'\n(%s) WARNING: input month character doesnt match any known 3-letter month name\n', mfilename);
end


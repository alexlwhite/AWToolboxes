%% function f = getNGramFreq(ngram, startYear, endYear, corpus, smoothing, caseInsenstive)
% This function queries Google Ngrams using the "webread" function to get
% the mean  frequency statistic for a single character string (ngram). It returns the frequencencies of each letter
% string averaged over a range of years. 
% 
% Inputs: 
% - ngram: a single character string you want to search for. For
%    example, "water fall" or "water"; 
% - startYear and endYear: the years over which we want to average 
% - corpus: a tag for an nGram cormpus. 'eng_2019' is a good bet. 
% - smoothing: how much smoothing to apply over years. Given that we're
%   averaging anyway, 0 is fine. 
% - caseInsensitive: character string "true" or "false". If true, then we
%   average over all case variations of the character string. 
% 
% Outputs: 
% - f: the mean ngram frequency, averaged over the date range provided. 
% 
% By Alex White, at Barnard College, 2025
function f = getNGramFreq(ngram, startYear, endYear, corpus, smoothing, caseInsenstive)

%case-insensitivity flag: accept either a logical or the strings 'true'/'false'
if ischar(caseInsenstive) || isstring(caseInsenstive)
    caseFlag = strcmpi(caseInsenstive, 'true');
else
    caseFlag = logical(caseInsenstive);
end
if caseFlag, caseStr = 'true'; else, caseStr = 'false'; end

%deal with spaces in the input ngram (URL-encode them):
ngram = strrep(ngram, ' ', '%20');

url = sprintf('https://books.google.com/ngrams/json?content=%s&year_start=%i&year_end=%i&corpus=%s&smoothing=%i&case_insensitive=%s', ngram, startYear, endYear, corpus, smoothing, caseStr);
data = webread(url);

if isempty(data) %empty data means nothing found? 
    f = NaN;
else
    if caseFlag && iscell(data)
    if strcmp(data{1}.type,'CASE_INSENSITIVE')
        data = data{1};
    else
        keyboard;
    end
end
   

f = mean(data.timeseries);

end
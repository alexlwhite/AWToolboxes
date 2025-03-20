%% function fs = getMultNGramsFreq(ngrams, startYear, endYear, corpus, smoothing, caseInsenstive)
% This function queries Google Ngrams using the "webread" function to get
% ngram frequency statistics. It returns the frequencencies of each letter
% string averaged over a range of years. 
% 
% Inputs: 
% - ngrams: a cell array of character strings you want to search for. For
%    example, {"water fall", "water"}; 
% - startYear and endYear: the years over which we want to average 
% - corpus: a tag for an nGram cormpus. 'eng_2019' is a good bet. 
% - smoothing: how much smoothing to apply over years. Given that we're
%   averaging anyway, 0 is fine. 
% - caseInsensitive: character string "true" or "false". If true, then we
%   average over all case variations of the character string. 
% 
% Outputs: 
% - fs: a vector of mean ngram frequencies, one for each character string
%   in "ngrams". 
%    If an item is not found, its frequency is set to NaN. 
% 
% By Alex White, at Barnard College, 2025


function fs = getMultNGramsFreq(ngrams, startYear, endYear, corpus, smoothing, caseInsenstive)

%deal with spaces in the input ngram:
spacer = '%20';
webOpts = weboptions('Timeout',25);

if ~iscell(ngrams)
    ngrams = {ngrams};
end

ngrams_sp = ngrams; %corrected with spacer inserted

query = '';
for ni=1:length(ngrams)
    if ni>1, query = [query ',']; end
    ngram = ngrams{ni};
    ispace = find(ngram==' ');
    if length(ispace)==1
        ngram = [ngram(1:(ispace-1)) spacer ngram((ispace+1):end)];
        ngrams_sp{ni} = ngram;
    elseif length(ispace)>1
        keyboard
    end

    query = [query ngram];
end


url = sprintf('https://books.google.com/ngrams/json?content=%s&year_start=%i&year_end=%i&corpus=%s=smoothing=%i&case_insensitive=%s', query, startYear, endYear, corpus, smoothing, caseInsenstive);
data = webread(url, webOpts);

fs = zeros(size(ngrams));

if isempty(data) %empty data means nothing found?
    fs = NaN;
else
    if iscell(data)
        names = cell(size(data));
        types = cell(size(data));

        for di=1:length(data)
            names{di} = data{di}.ngram;
            types{di} = data{di}.type;
        end


        for gi=1:length(ngrams)
            if caseInsenstive %find the sum of all possible versions with different letter cases
                thisOne = find(strcmp(names, [ngrams{gi} ' (All)']) & strcmp(types, 'CASE_INSENSITIVE'));
            else
                thisOne =  find(strcmp(names, ngrams{gi}));
            end
            if length(thisOne)==1
                fs(gi) = mean(data{thisOne}.timeseries);

            elseif isempty(thisOne)
                fs(gi)=NaN;
                fprintf('\nNGram %s not found\n', ngrams{gi});
            else
                keyboard
            end


        end



    else %for not case-insensitive requests, it comes back as a multiple-element struct
        nFound = numel(data);
        names = cell(1, nFound);
        for di=1:nFound
            names{di} = data(di).ngram;
        end
        for gi = 1:length(ngrams)
            thisOne =  find(strcmp(names, ngrams{gi}));
            if length(thisOne)==1
                fs(gi) = mean(data(thisOne).timeseries);

            elseif isempty(thisOne)
                fs(gi)=NaN;
                fprintf('\nNGram %s not found\n', ngrams{gi});
            else
                keyboard
            end

        end
    end
end

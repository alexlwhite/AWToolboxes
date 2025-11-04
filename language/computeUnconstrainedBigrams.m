function bigrams = computeUnconstrainedBigrams(word)
%bigrams are always put in lower case 

n = length(word);
bigrams = cell(1, n-1);
for i = 1:(n-1)
    bigrams{i} = lower(word([i, i+1]));
end

bigrams = unique(bigrams);


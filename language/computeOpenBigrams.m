function bigrams = computeOpenBigrams(word)
%an open bigram is a pair of letters that occur in the word, in a
%particular left-right order, regardless of any intervening letters. 
%these bigrams are stored in lower case always. 

n = length(word);
bigrams = {};
c = 0;
for i = 1:(n-1)
    for j=(i+1):n
        c = c+1;
        bigrams{c} = lower(word([i, j]));
    end
end

bigrams = unique(bigrams);


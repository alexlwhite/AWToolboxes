function bigrams = computeConstrainedBigrams(word)
%a constrained bigram is labeled as 'xy_i_j', where xy is the pair of
%letters (always put into lower case), i is the position of the 1st letter
%in the word, and j is the length of the word. 

n = length(word);
bigrams = cell(1,n-1);
for i = 1:(n-1)
    bigrams{i} = sprintf('%s_%i_%i',lower(word([i, i+1])), i, n);
end

bigrams = unique(bigrams);


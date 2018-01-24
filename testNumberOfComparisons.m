x = {'a','b','c','d','e','f','g','h','j'};

n = numel(x); 

pairs = {};
for i=1:n
    others = setdiff(1:n,i); 
    newpairs = cell(1,length(others));
    for j=1:length(others)
        newpairs{j} = sort([x{i} x{others(j)}]);
    end
    pairs = unique(cat(2,pairs,newpairs));
end

nPairs = numel(pairs)

predictedNPairs = (n-1)*n/2

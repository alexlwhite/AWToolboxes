function m = medianFromCounts(values, counts)

if length(values)~=length(counts)
    error('need the same number of values as counts');
end

xs = [];
for ii = 1:length(values)
    xs = [xs ones(1, counts(ii))*values(ii)];
end

m = median(xs);
function c = cosineSimilarity(x, y)


if ~isvector(x)|| ~isvector(y)
    error('x and y have to be vectors!')
end

if length(x)~=length(y)
    error('x and y have to be same length!')
end

if size(x,1)<size(x,2)
    x = x';
end
if size(y,1)<size(y,2)
    y = y';
end

c = dot(x, y)/(norm(x)*norm(y));
%y = matrixClone(x,a)
%
%by Alex White, 2013
%
%Takes a matrix x and returns a matrix y in which each element of x is
%repeated sequentially a times, along the horizontal (columns) dimension

%eg, x=[1 2 3; 
%       4 5 6];
%    a=3; 
%    y=[1 1 1 2 2 2 3 3 3; 
%       4 4 4 5 5 5 6 6 6]
%
%
% size(y) = [size(x,1) size(x,2)*a];


function y = matrixClone(x,a) 

len=size(x,2); 

y = zeros(size(x,1),len*a); 

for ri=1:size(x,1)
    for ai=1:a
        y(ri,ai:a:(ai+a*(len-1)))=x(ri,:);
    end
end
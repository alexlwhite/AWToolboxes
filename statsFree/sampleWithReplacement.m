%function y = sampleWithReplacement(n,k) 
%Inputs: 
% - n: either a vector containting the population to sample from, or a
% single interger, in which case we'll sample from the population 1:n. 
% - k: the number of samples to get 
%
% Outputs: a vector of length k 


function y = sampleWithReplacement(n,k) 

if numel(n)==1 && round(n) == n && n >= 0
    inputPopulation = false;
else
    inputPopulation = true;
    population = n;
    n = numel(population);
    if ~isvector(population)
        error(message('(sampleWithReplcement: input population n must be a vector)'));
    end
end

y = randi(n,k,1);
 
if inputPopulation
    y = population(y);
end
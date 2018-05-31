% function s1 = concatVectorsInStruct(s1,s2) 
% Alex White
% 
% concatenates vectors in two structures 
% 
% Inputs: 
% - s1: a structure with existing vectors as fields. 
% - s2: a similar structure with the fields with the same names, that are
% also vectors. 
% 
% Outputs
% - s1: the starting s1, with the data in each field of s2 concatenated on
% to each vector field of the original s1
%
% Note: if s1 lacks a field that s2 has, that field will be created in s1. 


function s1 = concatVectorsInStruct(s1, s2)

fs = fieldnames(s2);


for vi = 1:numel(fs)
    vn = fs{vi};
    if ~isfield(s1, vn)
        eval(sprintf('s1.%s = [];',vn));
    end
    eval(sprintf('s1.%s = [s1.%s s2.%s];',vn, vn, vn));
    
end

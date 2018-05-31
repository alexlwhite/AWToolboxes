function vns = fullFieldnames(s, vns) 

%By Alex White, April 2011
%Returns a list of all fieldnames in structure s. Recursively compiles
%them, so that if any field is itself a structure its fields will also be
%added to the list 

if nargin==1
    vns={}; 
end

pvs=fieldnames(s); 
for vi=1:numel(pvs)
    vn=pvs{vi}; 
    eval(sprintf('strct=isstruct(s.%s);', vn)); 
    if strct
        eval(sprintf('newvns=fullFieldnames(s.%s);', vn)); 
        for nvi=1:numel(newvns) 
            vns=cat(1,vns, sprintf('%s.%s', vn, newvns{nvi}));
        end
    else
        vns=cat(1,vns,vn);
    end
end

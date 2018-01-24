%function sems = withinSubjectSEMs(ds) 
%Compute within-subject SEMs following Morey (2008). 
%Subtract out each subject's global mean, then compute SEMs across subjects, 
%and correct with scale factor. 
%
%Input: ds: data matrix. Put separate conditions in columns and subjects in
%rows. 
%
%Output: sems: within-subject SEMs. Same size as ds, up to the last
%dimension. 

function sems = withinSubjectSEMs(ds) 

if ndims(ds)>2
    error('\nInput ds should have only 2 dimensions: conditions in columns and subjects in rows.\n');
end

nsubj = size(ds,1); %number of subjects 
ncond = size(ds,2); %number of conditions

%compute global mean for each subject 
globalMeans = mean(ds,2); 
dsNorm = ds - repmat(globalMeans,[1 ncond]); 

%scale factor to correct SEMs
scaleF = (ncond/(ncond-1));

%SEMs
sems = scaleF*std(dsNorm,0,1)/sqrt(nsubj);



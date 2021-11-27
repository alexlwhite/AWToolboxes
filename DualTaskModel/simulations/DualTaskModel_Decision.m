%% function Rs = DualTaskModel_Decision(Es, Cs)
% Generates responses in dual-task paradigm, by making a decision about
% sensory evidence E relative to criteria Cs. 
% 
% Inputs: 
% - Es, a 1xN vector of sensory evidence about N stimuli 
% - Cs, a 1xM vector of decision criteria, to produce a response from
% amount M-1 alternatives. Best if Cs(1)=-inf and Cs(2) = inf. 
% 
% Outputs: 
% - Rs, a 1xN vector of responses, as integers between 1 and (M-1). 
%
% by Alex White, 2017

function Rs = DualTaskModel_Decision(Es, Cs)

n = length(Es); 
Rs = NaN(1,n);

for i=1:n
    ds = Es(i)-Cs;  %compare E to criteria Cs
    ds(ds==0) = randn(1,1)*0.0001; %don't allow E to be exactly equal to a criterion 
    Rs(i) = find(diff(sign(ds))==-2); %find which criteria E is between
end

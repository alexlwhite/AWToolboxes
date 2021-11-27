%% function Cs = generateBalancedCriteria(Ms,SDs,nRespLevs)
% Given the set of possible means (Ms) and SDs that will generate sensory
% evidence in a simulated experiment, pick decision criteria Cs that will
% make unbiased responses when the simulated observer picks on each trial
% between nRespLevs possible responses. 
% 
% Inputs: 
% - Ms: set of means of Guassian distributions from which sensory evidence
% will be drawn in different conditions (we assume that all those
% conditions are equally likely)
% - SDs: SDs of those distributions 
% - nRespLevs: the number of response alternatives the simulated observer
% has (e.g., 4 confidence ratings). 
% 
% Outputs: 
% - Cs: a vector of length nRespLevs+1, equal to the criteria that separate
% ranges into which sensory evidence must fall to generate a particular
% response. Cs(1)=-inf; Cs(end)=inf. 
% 
% by Alex White, 2017

function Cs = generateBalancedCriteria(Ms,SDs,nRespLevs)

simRs = [];
for mi=1:length(Ms)
    for si=1:length(SDs)
        simRs = [simRs randn(1,10000)*SDs(si)+Ms(mi)];
    end
end

quants = 0:(1/nRespLevs):1;
Cs = quantile(simRs,quants);
%but set 1st C to negative infinity, and last to positive infinity, so that
%there are really only nRespLevs possible even if by chance we get a value
%of E more extreme than simulated here 
Cs(1)=-inf; Cs(end)=inf;
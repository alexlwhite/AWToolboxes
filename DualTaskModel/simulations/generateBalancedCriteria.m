%% function Cs = generateBalancedCriteria(Ms,SDs,nRespLevs)
% Given the set of possible means (Ms) and SDs that will generate sensory
% evidence in a simulated experiment, pick decision criteria Cs that will
% make unbiased responses when the simulated observer picks on each trial
% between nRespLevs possible responses. 
% 
% Inputs: 
% - Ms: set of means of Guassian distributions from which sensory evidence
% will be drawn in different conditions (we assume that all those
% conditions are equally likely). Can be a vector or matrix. 
% - SDs: SDs of those distributions. If SDs has the same number of elements
%   as Ms, then they are paired element-by-element (Ms(i) goes with SDs(i)).
%   Otherwise, every mean is combined with every SD. 
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
if numel(Ms)==numel(SDs)
    %means and SDs are paired element-by-element (use (:) so that ALL elements
    %of matrix inputs are included, e.g. both target-present and target-absent columns)
    for mi=1:numel(Ms)
        simRs = [simRs randn(1,10000)*SDs(mi)+Ms(mi)];
    end
else
    %otherwise, combine every mean with every SD
    for mi=1:numel(Ms)
        for si=1:numel(SDs)
            simRs = [simRs randn(1,10000)*SDs(si)+Ms(mi)];
        end
    end
end

quants = 0:(1/nRespLevs):1;
Cs = quantile(simRs,quants);
%but set 1st C to negative infinity, and last to positive infinity, so that
%there are really only nRespLevs possible even if by chance we get a value
%of E more extreme than simulated here 
Cs(1)=-inf; Cs(end)=inf;
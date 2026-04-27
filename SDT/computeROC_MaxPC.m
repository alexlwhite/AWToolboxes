% function [maxPC, maxI] = computeROC_MaxPC(HR, FR, propPresent) 
% Alex White, April 2015
% 
% This function computes the maximum proportion correct across points on an
% ROC curve defined by a set of (hit rate, false alarm rate) pairs. This
% was proposed by John Palmer to be a less biased measure of sensitivity
% from an ROC analysis, especialy in the context of dual-task performance
% and geneating predictions from the all-or-non-serial model. That model
% assumes that 1 stimulus is processed on each trial and the other is not
% at all. Different strategies for making a guess response about the
% stimulus that was not processed can have odd effects on the *area under*
% the ROC curve (Ag). But Palmer believes that this new measure, the
% maximum proportion correct across the ROC curve, is less biased and makes
% more straightforward predictions on an AOC. 
% 
% Given a set of points defined by the hit rates (hr) and false alarm rates (fr), or [hr, fr], for each one we calculate the
% total proportion correct by averaging the hit rate and the correct reject
% rate (1-fr). We can either assume that target-present and target-absent
% trials were equally likely, or adjust it according to the input
% propPresent that says the proportion of target-present trials for each
% point. 
% 
% Inputs: 
% - HR: a vector of length C, containing the hit rates computed at C
% different criteria, from conservative to liberal. HR should therefore be
% monotonically increasing. 
% - FR: a vector (also of length C) containing the corresponding false alarm rates. 
% - propPresent: a single number, the proportion of target-present trials
% in the experiment. Used to calcalate the proportion correct (pc) at each
% point along the ROC.
% 
% Outputs: 
% - maxPC: the maximum proportion correct across the C [hr, fr] pairs. For
% each pair i, pc(i) = propPresent*hr(i) + (1-propPresent)*(1-fr(i)). 
% - maxI: the index of the pair with the highest propoction correct. If
% multiple pairs have that same pc value, then maxI is the average of them.
% 


function [maxPC, maxI] = computeROC_MaxPC(hr, fr, propPresent)

if nargin<3
    propPresent = 0.5; 
end

if propPresent<0 || propPresent>1
    error('Prop present must be betwee 0 and 1\n');
end

pcs = propPresent*hr + (1-propPresent)*(1-fr);

maxPC = max(pcs);

maxI = mean(find(pcs==maxPC));


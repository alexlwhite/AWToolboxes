%% function Es = DualTaskModel_Encoding(Ms, Ss) 
% Simulate encoding stage of a single trial in dual-task paradigm. Simply
% draws random variables Es, which are the evidence for target presence. 
% 
% Input: 
% - Ms: means of the hypothetical Guassian distributions of E 
% - Ss: standard deviations of the hypothetical gaussian distributions of E
% 
% Output: 
% - Es: a matrix of same size as Ms, with values drawn from Gaussian
% distributions with means = Ms and SDs = Ss. Correspond to perceptual
% evidence for target presence. 
% 
% by Alex White, 2017

function Es = DualTaskModel_Encoding(Ms, Ss) 

Es = randn(size(Ms)).*Ss + Ms; 


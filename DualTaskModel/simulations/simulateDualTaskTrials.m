%% function Rs = simulateDualTaskTrials(Ms,SDs,Cs,switchSides)
% Simulates a block of trials, generating sensory evidence and choosing
% psychophysical responses. Inputs are the means and SDs of distributions
% of sensory evidence for 2 stimuli on each trial, as well as the criteria
% that those evidence values are compared against to generate a response.
% This function also implements "selection errors", in which the response
% to one side is based incorrectly on evidence for the other side.  
% 
% Inputs: 
% - Ms: a TxS matrix of the means of Gaussian distributions from which
% sensory evidence Es should be drawn, for each stimulus on each trial. T
% is the number of trials, S is the number of stimuli per trial. 
% - SDs: a TxS matrix of the SDs of those Gaussians. 
% - Cs: the decision criteria against which sensory evidence is judged to
% produce the responses recorded in Rs. 
% - switchSides: a Tx1 vector that determines whether the wrong side is
% responded to. If switchSides(t)==0, nothing is changed and all proceeds
% as it should. If switchSides(t)==1, then Es(2) is set to Es(1), meaning
% that the response for side 2 is erroneously based on evidence from side 1
% (but side 1 is responded to appropriately. Vice versa is switchSides(t)==2.  
%
% Outputs: 
% - Rs: a TxS matrix of responses the simulated observer made on each
% trial for each stimulus. 
%
% by Alex White, 2017

function Rs = simulateDualTaskTrials(Ms,SDs,Cs,switchSides)

%1. Generate internal responses, which are evidence for target presence: 
Es = DualTaskModel_Encoding(Ms,SDs);

%2. Sometimes, the wrong side is processed: 
Es(switchSides==2,1)=Es(switchSides==2,2);
Es(switchSides==1,2)=Es(switchSides==1,1);

%3. Make a decision and response about each stimulus, by comparing Es to criteria Cs. 
nT = size(Ms,1); %number of simulated trials 
Rs = NaN(size(SDs)); %subject's responses
for ti=1:nT
    Rs(ti,:) = DualTaskModel_Decision(Es(ti,:),Cs);
end
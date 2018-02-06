%% function [h,t] = makeHRF(n,k,d,maxT,TR,doPlot)
% Return a hemodynamic response function for MRI analysis, using G.
% Boynton's Gamma function and specified parameters 
%
% Inputs: 
% - n: phase delay (seconds). Increasing n delays the rise of the HRF.
% - k: time constant, a.k.a tau (seconds). Increasing k stretches out the HRF.
% - d: Pure delay pefore HRF begins to rise (seconds). a.k.a. delta. 
% - maxT: maximum number of seconds to include in the HRF. 
% - TR: length of the TR used in the MRI experiment being analyzed
%   (seconds). Each time point in the HRF increments by 1 TR seconds. 
% - doPlot: whether to plot the HRF (which happend in a new figure window)
% 
% Outputs: 
% - h: a Mx1 vector, describing the HRF at M timepoints (starting at 0). 
% - t: a Mx1 vector of timepoints, in units specified by 'units'

function [h,t] = makeHRF(n,k,d,maxT,TR,doPlot)

t = 0:TR:(maxT-TR);
t = t';

h = Gamma_gmb(n,k,t-d); 

if doPlot
    figure; 
    plot(t,h,'.-','Color',[0.4 0.4 0.8],'LineWidth',1,'MarkerSize',12); 
    xlabel('time (s)'); 
    ylabel('hemodynamic response'); 
end

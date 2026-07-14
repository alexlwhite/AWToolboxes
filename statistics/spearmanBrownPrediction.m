%% function rNew = spearmanBrownPrediction(rho, m)
% This function applies the "Spearman-Brown prediction formula" to
% "disattentuate" a correlation coeffecient. 
% It can be used to 'correct' a split-half reliability correlation
% coefficient to estimate would the correlation would be if you had
% actually run the full test twice. 
% https://en.wikipedia.org/wiki/Spearman%E2%80%93Brown_prediction_formula
% 
% Inputs: 
% - r: the original correlation coefficient(s). Can be a vector or matrix. 
% - m: the factor by which you are imagining you are scaling the amount of
%  data. If rho is based on splitting your data in half and correlating half
%  1 vs half 2, and you want to know what the reliability would be if you
%  had twice as much data, m=2. Defaults to 2. 
% 
% output: 
% - rNew: the disattentuated correlation coefficent:  m*r/(1+(m-1)*r);

function rNew = spearmanBrownPrediction(r, m)

if nargin<2
    m = 2; 
end
 
rNew = m*r./(1+(m-1)*r);
              
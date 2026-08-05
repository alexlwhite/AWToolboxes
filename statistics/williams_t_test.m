%% function [t_stat, p_val] = williams_t_test(rAB, rAC, rBC, N)
% Written with Gemini's help, August 5 2026, to test the null hypothesis that
% two dependent overlapping correlations are the same. For example, if you
% took 3 measures all from the same subjects, A, B and C. Now you want to
% test if A correlates with B as well as A correlatse with C. In other words, testing whether rAB=rAC. 
% You also need to input the correlation between B anc C, rBC. 
% 
% This is Williams' t-test, as advocated in James Steiger's 1980 paper:
% Steiger, J. H. (1980). Tests for comparing elements of a correlation matrix. Psychological bulletin, 87(2), 245.
% See his equation 7. 
% 
%Inputs: 
%- rAB: correlation between A and B
%- rAC: correlation between A and C
%- rBC: correlation between B and C
%- N: sample size (number of participants)
%
% Outputs: 
% -t-value: for the null hypothesis that rAB=rAC. 
% -pvalue: for the t-value (df=N-3). 

function [t_stat, p_val] = williams_t_test(rAB, rAC, rBC, N)
   
    R = [1,   rAB, rAC; 
         rAB, 1,   rBC; 
         rAC, rBC, 1];
    
    detR = det(R);
    r_bar = (rAB + rAC) / 2;
    
    num = (rAB - rAC) * sqrt((N - 1) * (1 + rBC));
    den = sqrt(2 * ((N - 1) / (N - 3)) * detR + r_bar^2 * (1 - rBC)^3);
    
    t_stat = num / den;
    df = N - 3;
    p_val = 2 * (1 - tcdf(abs(t_stat), df)); % Two-tailed p-value
end
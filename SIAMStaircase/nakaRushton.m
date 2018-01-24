function dPrime = nakaRushton(params, contrast) 

Rmax = params(1); 
C50  = params(2); 
n    = params(3); 


dPrime = (Rmax * (contrast.^n)) ./ ((contrast.^n) + C50.^n); 

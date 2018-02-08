function y=Gamma_gmb(n,k,t)
% GAMMA
%	y=Gamma(n,k,t)
%	returns a gamma function on vector t
%	y=(t/k).^(n-1).*exp(-t/k)/(k*factorial(n-1));
%	which is the result of an n stage leaky integrator.

%	6/27/95 Written by G.M. Boynton at Stanford University
%   4/19/11 Simplified it for Psychology 448/538 at U.W.
%
y = (t/k).^(n-1).*exp(-t/k)/(k*factorial(n-1));
y(t<0) = 0;
x=-5:0.01:5; 

mu = -2; 
var = 0.5; 
sigma = sqrt(var);

ncp = NormalCumulative(x,mu,var);
cdfp = normcdf(x,mu,sigma);
figure; plot(x,ncp,'g-','LineWidth',3); 
hold on; plot(x,cdfp,'r-')
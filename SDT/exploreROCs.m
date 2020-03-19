d = 2; 

muNS = 0; 
muS  = muNS+d;

c = -10:0.01:10;

fr = 1-normcdf(c,muNS,1);
cr = 1-fr;
hr = 1-normcdf(c,muS,1);

figure;
subplot(1,2,1); hold on; 
plot(hr,hr,'k--'); 
plot(fr,hr,'b-');
axis square;
xlabel('False alarm rate'); ylabel('Hit rate'); 

subplot(1,2,2); hold on; 
plot(hr,1-hr,'k--'); 
plot(cr,hr,'b-');
axis square;
xlabel('Correct reject rate'); ylabel('Hit rate'); 


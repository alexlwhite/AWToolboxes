dms = [1 4 5 2 1 ];
dss = [1 2 3 3 2];
nCond = length(dms);
nts = 150;
ds = zeros(nts, nCond);

ms = zeros(nCond, 1);
cis = zeros(nCond,2);

for ci=1:nCond
    ds(:,ci) = randn(nts,1)*dss(ci)+dms(ci);
    ms(ci) = mean(ds(:,ci));
    cis(ci, :) = boyntonBootstrap(@mean, ds(:,ci), 1000, 68);
end

MM = mean(ms)
MCI = mean(cis,1);

M = mean(ds(:))
CI = boyntonBootstrap(@mean, ds(:), 1000, 68);

CI3 = boyntonBootstrap(@mean, ms, 1000, 68)

figure; hold on;
plot(1, MM, 'bo');
plot([1 1], MCI, 'b-');
plot(2, M, 'ro');
plot([2 2], CI, 'r-');


xlim([0 3]); 
set(gca,'xtick',[1 2], 'XTickLabel',{'FirstByCond','NotByCond'});

diff(CI)/diff(MCI)
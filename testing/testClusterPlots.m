%TEST clusterPlots
clear; close all;

% Generate sample data for clustering
N = 100;
data = num2cell(rand(4, 2, N),3); % 100 points in 2D
s1 = size(data,1); s2 = size(data,2);

opt.connectLev1IndivPts = false;
opt.connectLev2IndivPts = false;

opt.meanSymbol = '-'; %'-' means a line
opt.meanLineWidth = 5;

opt.xTickLabs = {'A','B','C','D'};
opt.yLab = 'madeup';
opt.symmetricErrorBar = false;
opt.errorBarCI = 95;
opt.legendLabs = {'XY','YZ'};

c1 = hsv2rgb([0.3 0.6 0.9]);
c2 = hsv2rgb([0.3 0.8 0.5]);
c3 = hsv2rgb([0.3 1.0 0.3]);
c4 = hsv2rgb([0.3 1.0 0.3]);

opt.fillColors = ones(s1,s2,3);

l2colrs = [c1; c2];
opt.edgeColors = repmat(reshape(l2colrs,[1 s2 3]), s1);




figure; hold on;

%plot with my old version 
[barCenters, opt, allX] = clusterPlot(data, opt);

%check differences between each subject's jitters
jitts = NaN(s1*s2, N);
cc = 0;
for i1 = 1:s1
    for i2=1:s2
        cc = cc+1;
        jitts(cc, :) = allX{i1, i2, :} - barCenters(i1, i2);
    end
end

jitSD = std(jitts, 0, 1);

maxSD1 = max(abs(jitSD))



figure; hold on;

c1 = hsv2rgb([0.7 0.6 0.8]);
c2 = hsv2rgb([0.7 0.8 0.5]);
c3 = hsv2rgb([0.7 1.0 0.3]);

opt.fillColors = ones(s1,s2,3);

l2colrs = [c1; c2];
opt.edgeColors = repmat(reshape(l2colrs,[1 s2 3]), s1);

[barCenters, opt, allX] = clusterPlot_Gemini4(data, opt);

%check differences between each subject's jitters
jitts2 = NaN(s1*s2, N);
cc = 0;
for i1 = 1:s1
    for i2=1:s2
        cc = cc+1;
        jitts2(cc, :) = allX{i1, i2, :} - barCenters(i1, i2);
    end
end

jitSD2 = std(jitts2, 0, 1);
maxSD2 = max(abs(jitSD2))

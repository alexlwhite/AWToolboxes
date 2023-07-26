%% Example script to use AWToolbox's barPlot_AW function 

%make up some data

%means: a 3x2 matrix
Ms = [0.6 0.7; 
      0.7 0.65; 
      0.8 0.9]; 

%standard errors: a 3x2 matrix
Es = randn(size(Ms))*0.1;

%set options for plot
opt.ylims = [0.33 1]; %y axis limits
opt.yticks = 0.4:0.2:1; %y axis tick marks
opt.xTickLabs = {'Left','Center','Right'}; %labels for x-axis ticks
opt.xLab = 'Word location'; %x-axis label
opt.yLab = 'Proportion correct'; %y-axis label
opt.legendLabs = {'Yes','No'}; %legend labels
opt.legendTitle = 'Has right VWFA';
opt.legendLoc = 'NorthWest'; %where the legend goes 

%set colors: for the two conditions that are in columns in the means
yesColr = [255,73,89]/255;
noColr = [0, 105, 146]/255;

%reshape colors to be in the right format (nRows x nCols x 3)

%edgeColors are the outlines of the bars
opt.edgeColors = NaN(size(Ms,1), size(Ms,2), 3);
opt.edgeColors(:, 1, :) = repmat(yesColr, 3, 1, 1); 
opt.edgeColors(:, 2, :) = repmat(noColr, 3, 1, 1); 

%fillColors are how the bas are filled
opt.fillColors = opt.edgeColors;

%error bar colors (here, all black)
opt.errorBarColors = zeros(size(opt.fillColors));

%now plot! 
figure; hold on;
barCtrs = barPlot_AW(Ms, Es, opt);
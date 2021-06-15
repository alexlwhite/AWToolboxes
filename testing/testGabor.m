width = 5;
height = 5; 
sdx = 1;
sdy = 0.75;
sf = 2;
angle = 30;
phase = 0;
contrast = 0.5;
xCenter = 0; 
yCenter = 0; 
bgLum = 0.5; 
deltaLum = 0.5; 
pixPerDeg = 100;

img = makeGabor(width, height, sdx, sdy, sf, angle, phase, contrast, xCenter, yCenter, bgLum, deltaLum, pixPerDeg);

figure; imshow(img)

angle = 90;

%% test gratings:
%close all;

contrast = 0.05;
bgLum = 0.6;
deltaLum = 0.4;
sf = 0.5;
g1 = makeGrating(width, NaN, sf, angle, phase, pixPerDeg);
g1 = g1*contrast*deltaLum + bgLum;


%Clip
if bgLum>1 && deltaLum>1
    white = 255;
else
    white = 1;
end
black = 0;

g1(g1>white) = white;
g1(g1<black) = black;

figure; subplot(2,1,2); 
plot(1:length(g1), g1);
xlim([1 length(g1)]);
ylim([black white]);
axis square;
axis off

g2 = makeGrating(width, height, sf, angle, phase, pixPerDeg);
g2 = contrast*g2*deltaLum + bgLum;
g2(g2>white) = white;
g2(g2<black) = black;


 subplot(2,1,1); 
 imshow(g2);

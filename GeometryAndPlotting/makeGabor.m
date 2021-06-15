%% img = makeGabor(width, height, sdx, sdy, sf, angle, phase, contrast, xCenter, yCenter, bgLum, deltaLum, pixPerDeg)
%      usage: img = makeGabor(width, height, sdx, sdy, sf, angle, phase, contrast, xCenter, yCenter, bgLum, deltaLum, pixPerDeg)
%         by: Alex White  
%       date: 1/11/21
%    purpose: create a Gabor patch.  
%    inputs: 
%             - width, height: size of the output image in degrees of
%             visual angle. If height is NaN then we just get a 1D
%             sinusoid. 
%             - sdx and sdy: standard deviations of Gaussian window in horizontal and vertical directions, in degrees of visual angle
%             - sf is spatial frequency in cycles/degree of visual angle 
%             - angle: grating angle in degrees with 0 being horizontal 
%             - phase: phase of grating in degrees. 
%             - contrast: between 0 and 1
%              - xcenter and ycenter are the coordinates of the Gaussian's  center 
%            - bgLum: brightness level of the background
%            - deltaLum: max amount by which brightness can vary from bgLum
%             - pixPerDeg is the number of pixels per degree of visual angles. 

function img = makeGabor(width, height, sdx, sdy, sf, angle, phase, contrast, xCenter, yCenter, bgLum, deltaLum, pixPerDeg)

%sinusoidal grating, -1<grat<1
grat = makeGrating(width, height, sf, angle, phase, pixPerDeg);

gauss = make2DGaussian(width,height,sdx,sdy,xCenter,yCenter,pixPerDeg);

g = contrast*gauss.*grat;

img = g*deltaLum + bgLum;


%Clip
if bgLum>1 && deltaLum>1
    white = 255;
else
    white = 1;
end
black = 0;

img(img>white) = white;
img(img<black) = black;




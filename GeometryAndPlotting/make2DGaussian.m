% g = make2DGaussian(width,height,sdx,sdy,xCenter,yCenter,pixPerDeg)
%
%      usage: make2DGaussian(width,height,sdx,sdy,xCenter,yCenter,pixPerDeg)
%         by: Alex White (adopted from mglMakeGraussian by Justin Gardner)
%       date: 1/11/21
%    purpose: make a 2D gaussian, useful for making gabors. 
%    inputs: 
%             - width, height: size of the output image in degrees of visual angle
%             - sdx and sdy: standard deviations of Gaussian in horizontal and vertical directions, in degrees of visual angle
%             - xcenter and ycenter are the coordinates of the Gaussian's  center 
%             - pixPerDeg is the number of pixels per degree of visual angles. 
%    output: g, a height*pixPerDeg x width*pixPerDeg matrix, in the range [0, 1]

function g = make2DGaussian(width,height,sdx,sdy,xCenter,yCenter,pixPerDeg)

% get size in pixels
widthPixels = round(width*pixPerDeg);
heightPixels = round(height*pixPerDeg);
widthPixels = widthPixels + mod(widthPixels+1,2);
heightPixels = heightPixels + mod(heightPixels+1,2);

% get a grid of x and y coordinates that has 
% the correct number of pixels
x = -width/2:width/(widthPixels-1):width/2;
y = -height/2:height/(heightPixels-1):height/2;
[xMesh,yMesh] = meshgrid(x,y);

% compute gaussian window
g = exp(-(((xMesh-xCenter).^2)/(2*(sdx^2))+((yMesh-yCenter).^2)/(2*(sdy^2))));
% clamp small values to 0 so that we fade completely to gray.
g(g(:)<0.01) = 0;
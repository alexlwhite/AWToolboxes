% g = makeGrating(width,height,sf,angle,phase,pixPerDeg)
%
%      usage: makeGrating(width,height,sf,angle,phase,pixPerDeg)
%         by: Alex White (adopted from mglMakeGrating by Justin Gardner)
%       date: 1/11/21
%    purpose: create a 2D grating.  
%    inputs: 
%             - width, height: size of the output image in degrees of
%             visual angle. If height is NaN then we just get a 1D
%             sinusoid. 
%             - sf is spatial frequency in cycles/degree of visual angle 
%             - angle: grating angle in degrees with 0 being horizontal 
%             - phase: phase of grating in degrees. 
%             - pixPerDeg is the number of pixels per degree of visual angles. 
% 
%    output: g, a height*pixPerDeg x width*pixPerDeg matrix, ranging from -1 to 1

function g = makeGrating(width, height, sf, angle, phase, pixPerDeg)


% make it so that angle of 0 is horizontal
angle = -90-angle; 


% get size in pixels
widthPixels = round(width*pixPerDeg);
heightPixels = round(height*pixPerDeg);
widthPixels = widthPixels + mod(widthPixels+1,2);
heightPixels = heightPixels + mod(heightPixels+1,2);

% calculate image parameters
phase = pi*phase/180;

% if height is nan, it means we should calculate a 1 dimensional grating
if isnan(height)
  % 1D grating (note we ignore orientation)
  x = -width/2:width/(widthPixels-1):width/2;
  g = cos(x*sf*2*pi+phase);
else
  % 2D grating
  % calculate orientation
  angle = pi*angle/180;
  a=cos(angle)*sf*2*pi;
  b=sin(angle)*sf*2*pi;

  % get a grid of x and y coordinates that has 
  % the correct number of pixels
  x = -width/2:width/(widthPixels-1):width/2;
  y = -height/2:height/(heightPixels-1):height/2;
  [xMesh,yMesh] = meshgrid(x,y);

  % compute grating
  g = cos(a*xMesh+b*yMesh+phase);
end
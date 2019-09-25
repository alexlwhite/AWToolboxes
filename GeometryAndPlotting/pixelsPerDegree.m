function [pixelsPerDegree] = pixelsPerDegree(distance, width, numPixels)
% [pixelsPerDegree] = pixelsPerDegree(distance, width, numPixels)
%
% First computes the size of on monitor pixel into degrees of visual angle with the formula:
%
% ang = 2 * arctan(pixSizeCm/(2*distance)) * (180/pi);
% 
% Then computes the number of pixels in 1 degree of visual angle: 
% pixelsPerDegree = 1/ang;
%
% Input:
%       distance        Distance from screen, cm
%       width           Width of screen, cm
%       numPixels       Screen resolution in pixels, width 
% 
%
% Output: 
%   pixelsPerDegree     number of pixels in 1 degree of visual angle 
%
% 
% Note:
% - Warning: Assumes isotropic (square) pixels

% By Alex White Adpoted from pix2Deg by G.M. Boynton, Zach Ernst & Kely Chang - 11/1/07
% 
%% Convert Pixels to Visual Angels

%size of one pixel in centimeters
pixSize = width / numPixels; % cm/pix

%size of one pixel in degrees of visual angle 
ang = 2*atan(pixSize/(2*distance))*(180/pi);

%pixels per 1 degree of visual angle
pixelsPerDegree = 1/ang;

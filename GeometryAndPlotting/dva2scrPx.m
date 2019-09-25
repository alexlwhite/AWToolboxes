% [scrXY] = dva2scrPx(scr, x,y)
% 
% given coordinates in degrees of visaul angle with respect to the 
% screen center, with positive numbers being up and to the right, 
% returns coordinates on the screen in terms of pixels from upper left

function [scrXY] = dva2scrPx(scr, x,y)


xpx=x*scr.ppd;
ypx=y*scr.ppd; 

scrXY=[scr.centerX+xpx scr.centerY-ypx];



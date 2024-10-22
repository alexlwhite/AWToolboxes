%% function ys = wrapAroundTilts(xs)
% When simulating the degrees of tilt of a line or Gabor patch, it gets
% complicated because after you rotate a line by 180 deg, youre back to
% where youre started. And when you go past 90 of tilt, you're equivalently
% at some tilt between -90 and 0. 
% This function takes care of that, taking any tilts in degrees from
% -infinity to infinity, and wrapping them back around, that is,
% constraining them to the range [-90 90]

function ys = wrapAroundTilts(xs)

%make these values wrap around the circle
ys = mod(xs, 180);

%now deal with the fact that 95deg is acutally the same as -85
ys(ys>90) = ys(ys>90)-180;

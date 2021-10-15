%% function bounds = ImageBoundsAW(img, bgColor) 
% Return tight bounding box within an image that's not background. 
% Based on PsychToolbox's ImageBounds, which seems to have starting rows
% and starting columns that include 1 background pixel. 
% 
% inputs: 
% - img: an RxCx3 or RxCx1 image. If the image is RGB (size(img,3)>1), then
% we average over the 3rd dimension before finding non-background pixels. 
% 
% - bgColor: a single number that states the *grayscale* value of the background color. This
% function finds the bounding box of the image that's not background. RGB
% colors of background are currently not implemented. But could be. 
% 
% Outputs: 
% bounds: [left top right bottom]  indices of the maximal extent of the non-background image. 
% bounds(1), the left edge, is the first column that has any image in it. 
% bounds(3), the right edge, is the last column that has any image in it. 
% bounds(2), the top edge, is the first row that has any image in it. 
% bounds(4), the bottom edge, is the last row that has any image in it. 

function bounds = ImageBoundsAW(img, bgColor) 

if nargin<2
    bgColor = 255;
end

% Find all non-background pixels:
if length(bgColor)==1
   if size(img,3)>1
      img=mean(img,3);
   end
   [rows,cols]=find(img(:,:)~=bgColor);
else
   error('Support for color images not yet implemented.');
end

% Compute their bounding rect and return it:
if isempty(rows) || isempty(cols) %image is empty
    bounds=[0 0 0 0];
else
    % 	newRect=[left,top,right,bottom];
    bounds=[min(cols), min(rows), max(cols), max(rows)];
end
return;

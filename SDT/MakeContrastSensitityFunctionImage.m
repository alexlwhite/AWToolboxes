% demoMonoCSF

% Draws a variant of a Campbell-Robson CSF chart using a Bits++.

% -----------------------------------------------------------------------------

% This demonstration illustrates the use of the high bit depth drawing modes

% available on the Bits++ to display (a variant of) the Campbell-Robson

% Contrast Sensitivity Function chart. This particular demo uses the high

% bit depth mono drawing mode. (14bitGreyScaleMode). Another similar

% demonstration is available that uses the high resolution colour mode

% (42bitColourMode).

%

% View the source code for this file for pointers on how to use the high bit

% depth mono drawing mode.



% Please be aware that in this source description, some of the CSF chart

% parameters have been modified to improve it's appearance on screen.

% Please do not rely on this file as a reference for how to produce a CSF

% chart.



% originally written by wp for Visage

% converted to Bits++ and OSX Psychtoolbox by ejw



clear; close all


% 
% % Define screen
% 
% whichScreen=max(Screen('Screens'));
% 
% 
% 
% % Find the color values which correspond to white and black.  Though on OS
% 
% % X we currently only support true color and thus, for scalar color
% 
% % arguments,
% 
% % black is always 0 and white 255, this rule is not true on other platforms will
% 
% % not remain true on OS X after we add other color depth modes.  
% 
% white=WhiteIndex(whichScreen);
% 
% black=BlackIndex(whichScreen);
% 
% gray=GrayIndex(whichScreen);

white = 255; 
black = 0; 
gray = 127.5;


% Open up the a window on the screen.

% This is done by calling the Screen function of the Psychophysics

% Toolbox. We need to open the onscreen window first, so that

% the offscreen storage allocated subsequently is properly

% matched to the onscreen frame buffer.

%[window,screenRect] = Screen('OpenWindow',whichScreen,128,[],32,2);



% THE FOLLOWING STEP IS IMPORTANT.

% make sure the graphics card LUT is set to a linear ramp

% (else the encoded data will not be recognised by Bits++).

% There is a bug with the underlying OpenGL function, hence the scaling 0 to 255/256.  

% This demo will not work using a default gamma table in the graphics card,

% or even if you set the gamma to 1.0, due to this bug.

% This is NOT a bug with Psychtoolbox!

%Screen('LoadNormalizedGammaTable',window,linspace(0,(255/256),256)'*ones(1,3));



% draw a gray background on front and back buffers
% 
% Screen('FillRect',window, gray);
% 
% Screen('Flip', window);
% 
% Screen('FillRect',window, gray);
% 


% =================================================================

% CODE NEEDED HERE !

% "linear_lut" should be replaced here with one giving the inverse

% characteristic of the monitor.

% =================================================================

% restore the Mono++ overlay LUT to a linear ramp

%linear_lut =  repmat(round(linspace(0, 2^16 -1, 256))', 1, 3);

%BitsPlusSetClut(window,linear_lut);



% find out how big the window is

%[POINTS_ACROSS, POINTS_DOWN]=Screen('WindowSize', window);

POINTS_ACROSS = 1280;
POINTS_DOWN = 960;

% ===========================

% DEFINE SPATIAL FREQUENCY GP

% ===========================

SF_LOW        = 1.5;

SF_HIGH       = 400;

SF_BASE       = (SF_HIGH/SF_LOW) .^ (1/POINTS_ACROSS); % base

SF_POW        = [1:POINTS_ACROSS];                     % power

SF_FORM       = SF_BASE .^ SF_POW;                     % functional form

SF_MUL        = (SF_LOW/SF_BASE);                      % scaling factor

SF_VALU       = SF_FORM * SF_MUL;                      % values



% ==================

% DEFINE CONTRAST GP

% ==================

CT_LOW      = 1/1024;

CT_HIGH     = 1;

CT_BASE     = (CT_HIGH/CT_LOW) .^ (1/POINTS_DOWN); % base

CT_POW      = [1:POINTS_DOWN];                     % power

CT_FORM     = CT_BASE .^ CT_POW;                   % functional form

CT_MUL      = (CT_LOW/CT_BASE);                    % scaling factor

CT_VALU     = CT_FORM .* CT_MUL;                   % values



% ============================

% FORM IMAGE FROM PROGRESSIONS

% ============================

[XGRID,YGRID] = meshgrid(sin(SF_VALU),CT_VALU); % CT & SF are orthogonal

IMAGE         = XGRID .* YGRID;                 % contrast modulates amplitude

IMAGE         = (IMAGE + 1) / 2;                % map values to between 0 and 1.



% ==================================================================

% CODE NEEDED HERE !

% "IMAGE" should be corrected here for the inverse characteristic of

% the monitor.

% ================================================================== 



% encode the image in mono++ format

%encoded_image = BitsPlusPackMonoImage(IMAGE*((2^16)-1));

  

% put the csf chart image up on the back buffer

figure; imshow(IMAGE); 


% Screen('PutImage', window, IMAGE);
% 
% 
% 
% % make the make buffer current during the next blanking period
% 
% Screen(window,'Flip');
% 
% 
% 
% fprintf('Displaying CSF chart, hold a key to exit \n');
% 
% 
% 
% while(~KbCheck)
% 
% end;    


% 
% % if the system only has one screen, set the LUT in Bits++ to a linear ramp
% 
% % if the system has two or more screens, then blank the screen.
% 
% if (whichScreen == 0)
% 
%     % =================================================================
% 
%     % CODE NEEDED HERE !
% 
%     % "linear_lut" should be replaced here with one giving the inverse
% 
%     % characteristic of the monitor.
% 
%     % =================================================================    
% 
%     % restore the Mono++ overlay LUT to a linear ramp
% 
%     %linear_lut =  repmat(round(linspace(0, 2^16 -1, 256))', 1, 3);
% 
%  %   BitsPlusSetClut(window,linear_lut);
% 
%     
% 
%     % draw a gray background on front and back buffers to clear out any old LUTs
% 
%     Screen('FillRect',window, gray);
% 
%     Screen('Flip', window);
% 
%     Screen('FillRect',window, gray);
% 
%     Screen('Flip', window);
% 
% 
% 
%     % Close the window.
% 
%     Screen('CloseAll');    
% 
% else
% 
%     % Blank the screen (zero the Mono++ overlay LUT),
% 
%     % - doesn't actually quite blank the screen.
% 
% %    BitsPlusSetClut(window,zeros(256,3));
% 
% 
% 
%     % draw a black background on front and back buffers to clear out any old LUTs
% 
%     Screen('FillRect',window, black);
% 
%     Screen('Flip', window);
% 
%     Screen('FillRect',window, black);
% 
%     Screen('Flip', window);
% 
% 
% 
%     Screen('CloseAll');
% 
% end
% 

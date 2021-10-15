function [letterStruct,alphabetBounds]=CreateLetterTextures_AW(condition,o,window)
% [letterStruct,alphabetBounds]=CreateLetterTextures(condition,o,window)
% Create textures, one per letter in o.alphabet plus o.borderLetter.
% Returns array "letterStruct" with one struct element per letter, plus
% bounding box "alphabetBounds" that will hold any letter. Called by
% CriticalSpacing.m. The font is o.targetFont.
%
% If o.readAlphabetFromDisk==0 then the font is rendered by Screen DrawText
% to create a texture for each desired letter. The font's TextSize is
% computed to yield the desired o.targetPix size in the direction specified
% by o.targetSizeIsHeight. However, if
% o.targetFontHeightOverNominalPtSize==nan then the TextSize is set equal
% to o.targetPix.
%
% If o.readAlphabetFromDisk==1 then we look for a folder inside
% CriticalSpacing/lib/alphabets/ whose name matches that of the desired
% font. We give a fatal error if it's not found. The folder is very simple,
% one image file per letter; the filename is the letter, URL encoded to
% allow you to specify symbols, including a space.
%
% The argument "condition" is used only for diagnostic printout.
%
% HISTORY 
% Feb 22, 2017 - Alex White (ALW) cloned this from Pelli's github
% Feb 22, 2017 - ALW added functionality to also return image of each word,
% as a matrix of pixel values. That's now in letterStruct(i).image. Also
% returns final pixel size of font, in letterStruct(i).sizePix;

if ~isfinite(o.targetHeightOverWidth)
   o.targetHeightOverWidth=1;
end
letters=[o.alphabet o.borderLetter];
for i=1:length(letters)
   letterStruct(i).letter=letters(i);
end
% if o.targetSizeIsHeight
%    canvasRect=[0 0 o.targetPix o.targetPix];
% else
   canvasRect=[0 0 o.targetPix o.targetPix]*o.targetHeightOverWidth;
% end
black=0;
white=o.bgColor; %was 255
if o.readAlphabetFromDisk

   % Read from disk into "savedAlphabet".
   alphabetsFolder=fullfile(fileparts(fileparts(mfilename('fullpath'))),'alphabets'); % CriticalSpacing/alphabets/
   if ~exist(alphabetsFolder,'dir')
      error('Folder missing: "%s"',alphabetsFolder);
   end
   folder=fullfile(alphabetsFolder,urlencoding(o.targetFont));
   if ~exist(folder,'dir')
      error('Folder missing: "%s". Target font "%s" has not been saved.',folder,o.targetFont);
   end
   d=dir(folder);
   ok=~[d.isdir];
   for i=1:length(ok)
      systemFile=streq(d(i).name(1),'.') && length(d(i).name)>1;
      systemFile=systemFile || streq(d(i).name,'Thumbs.db'); % ignore Windows cache
      ok(i)=ok(i) && ~systemFile;
   end
   d=d(ok);
   if length(d)<length(o.alphabet)
      error('Sorry. Saved %s alphabet has only %d letters, and you requested %d letters.',o.targetFont,length(d),length(o.alphabet));
   end
   savedAlphabet.letters=[];
   savedAlphabet.images={};
   savedAlphabet.rect=[];
   for i=1:length(d)
      filename=fullfile(folder,d(i).name);
      try
         savedAlphabet.images{i}=imread(filename);
      catch
         sca;
         error('Cannot read image file "%s".',filename);
         psychrethrow(psychlasterror);
      end
      if isempty(savedAlphabet.images{i})
         error('Cannot read image file "%s".',filename);
      end
      [~,name]=fileparts(urldecoding(d(i).name));
      if length(name)~=1
         error('Saved "%s" alphabet letter image file "%s" must have a one-character filename after urldecoding.',o.targetFont,name);
      end
      savedAlphabet.letters(i)=name;
      savedAlphabet.bounds{i}=ImageBounds(savedAlphabet.images{i},255);
      savedAlphabet.imageBounds{i}=RectOfMatrix(savedAlphabet.images{i});
      if o.showLineOfLetters
         fprintf('%d: CreateLetterTextures "%c" image(%d) width %d, ',condition,savedAlphabet.letters(i),i,RectWidth(savedAlphabet.bounds{i}));
         fprintf('bounds %d %d %d %d, image %d %d %d %d.\n',savedAlphabet.bounds{i},savedAlphabet.imageBounds{i});
      end
      if isempty(savedAlphabet.rect)
         savedAlphabet.rect=savedAlphabet.bounds{i};
      else
         a=savedAlphabet.rect;
         savedAlphabet.rect=UnionRect(savedAlphabet.rect,savedAlphabet.bounds{i});
         b=savedAlphabet.rect;
         if ~all(a==b)
           warning('Found a change in alphabet size!'); % put a breakpoint here
           % keyboard
         end
      end
   end
   alphabetBounds=savedAlphabet.rect;

   % Create textures, one per letter.
   for i=1:length(letters)
       which=strfind([savedAlphabet.letters],letters(i));
      if length(which)~=1
         error('Letter %c is not in saved "%s" alphabet "%s".',letters(i),o.targetFont,savedAlphabet.letters);
      end
      assert(length(which)==1);
      r=savedAlphabet.rect;
      letterImage=savedAlphabet.images{which}(r(2)+1:r(4),r(1)+1:r(3));
      if o.contrast==1
         letterStruct(i).texture=Screen('MakeTexture',window,letterImage);
      else
         letterStruct(i).texture=Screen('MakeTexture',window,uint8(white+(double(letterImage)-white)*o.contrast));
      end
      letterStruct(i).rect=Screen('Rect',letterStruct(i).texture);
      % Screen DrawTexture will later scale and stretch, as needed.
   end

else % if o.readAlphabetFromDisk
   % Draw font and get bounds.
   scratchWindow=Screen('OpenOffscreenWindow',window,[],canvasRect*4,8,0);
   if ~isempty(o.targetFontNumber)
      Screen('TextFont',scratchWindow,o.targetFontNumber);
      [~,number]=Screen('TextFont',scratchWindow);
      assert(number==o.targetFontNumber);
   else
      Screen('TextFont',scratchWindow,o.targetFont);
      font=Screen('TextFont',scratchWindow);
      assert(streq(font,o.targetFont));
   end
   if o.targetSizeIsHeight
      sizePix=round(o.targetPix/o.targetFontHeightOverNominalPtSize);
   else
      sizePix=round(o.targetPix*o.targetHeightOverWidth/o.targetFontHeightOverNominalPtSize);
   end
   if ~isfinite(sizePix)
      sizePix=o.targetPix;
   end
   Screen('TextSize',scratchWindow,sizePix);
   for i=1:length(letters)
      lettersInCells{i}=letters(i);
      bounds=TextBoundsNew(scratchWindow,letters(i),1);
      if o.showLineOfLetters
         b=Screen('TextBounds',scratchWindow, letters(i));
         fprintf('%d: %s "%c" textSize %d, TextBounds [%d %d %d %d] width x height %d x %d, Screen TextBounds %.0f x %.0f\n', ...
            condition,o.targetFont,letters(i),sizePix,round(bounds),RectWidth(bounds),RectHeight(bounds),RectWidth(b),RectHeight(b));
      end
      letterStruct(i).bounds=bounds;
      if i==1
         alphabetBounds=bounds;
      else
         alphabetBounds=UnionRect(alphabetBounds,bounds);
      end
   end
   bounds=alphabetBounds;
   if o.printSizeAndSpacing
      fprintf('%d: sizePix %d, first letter "%c", height %d, width %d.\n',condition,sizePix,letters(1),RectHeight(letterStruct(1).bounds),RectWidth(letterStruct(1).bounds));
   end
   assert(RectHeight(bounds)>0);
   for i=1:length(letters)
      letterStruct(i).width=RectWidth(letterStruct(i).bounds);
      desiredBounds=CenterRect(letterStruct(i).bounds,bounds);
      letterStruct(i).dx=desiredBounds(1)-letterStruct(i).bounds(1);
   end
   Screen('Close',scratchWindow);

   % Create texture for each letter
   canvasRect=bounds;
   canvasRect=OffsetRect(canvasRect,-canvasRect(1),-canvasRect(2));
   if o.printSizeAndSpacing
      fprintf('%d: textSize %.0f, "%s" height %.0f, width %.0f\n',condition,sizePix,letters,RectHeight(bounds),RectWidth(bounds));
   end
   for i=1:length(letters)
      [letterStruct(i).texture,letterStruct(i).rect]=Screen('OpenOffscreenWindow',window,[],canvasRect,8,0);
      if ~isempty(o.targetFontNumber)
         Screen('TextFont',letterStruct(i).texture,o.targetFontNumber);
         [font,number]=Screen('TextFont',letterStruct(i).texture);
         assert(number==o.targetFontNumber);
      else
         Screen('TextFont',letterStruct(i).texture,o.targetFont);
         font=Screen('TextFont',letterStruct(i).texture);
         assert(streq(font,o.targetFont));
      end
      Screen('TextSize',letterStruct(i).texture,sizePix);
      Screen('FillRect',letterStruct(i).texture,white);
      if o.contrast==1
         Screen('DrawText',letterStruct(i).texture,letters(i),-bounds(1)+letterStruct(i).dx,-bounds(2),black,white,1);
      else
         Screen('DrawText',letterStruct(i).texture,letters(i),-bounds(1)+letterStruct(i).dx,-bounds(2),uint8(white+(double(black)-white)*o.contrast),white,1);
      end
      
      %CAPTURE IMAGE!
      letterStruct(i).image=double(Screen('GetImage',letterStruct(i).texture,[],[],0,3));
      letterStruct(i).sizePix = sizePix;
      
   end
end % if o.readAlphabetFromDisk
end

function u = urlencoding(s)
  u = '';

  for k = 1:length(s),
    if ~isempty(regexp(s(k), '[a-zA-Z0-9]', 'once'))
      u(end+1) = s(k);
    else
      u=[u,'%',dec2hex(s(k)+0)];
    end;
  end
end

function u = urldecoding(s)
  u = '';
  k = 1;
  while k<=length(s)
    if s(k) == '%' && k+2 <= length(s)
      u = sprintf('%s%c', u, char(hex2dec(s((k+1):(k+2)))));
      k = k + 3;
    else
      u = sprintf('%s%c', u, s(k));
      k = k + 1;
    end
  end
end

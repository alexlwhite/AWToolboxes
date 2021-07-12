key = '.';

KbName('UnifyKeyNames');

keyNum = KbName(key); 

fprintf(1,'\nWaiting for key %s\n', key);
keyPress = 0;
while ~keyPress
    [keyPress] = checkTarPress(keyNum);
end
fprintf(1,'\nKey %s worked!\n', key);


%function esca
% executes "sca" for "screen('Close all') *and* then tries to stop eyelink
% recording and shutdown the eyelink connection 
function esca

sca;

if Eyelink('IsConnected')
    Eyelink('stoprecording');
    Eyelink('closefile');
    Eyelink('shutdown');
end
function [msac, monol, monor] = mergeBinocSaccades(sacl, sacr, mergeInt)
% mergeBinocSaccades Determination of binocular saccades from tables
%
% INPUTS:
%   sacl: Table from detectSaccades (left eye)
%   sacr: Table from detectSaccades (right eye)
%   mergeInt: the minimum separation (in samples) between an offset and the next
%             onset for the two saccades to be counted as separate. Otherwise they are
%             merged.
%
% OUTPUTS:
%   msac: Table of binocular saccades (20+ columns: averaged variables, then right eye variables, then left)
%   monol: Table of monocular left eye saccades
%   monor: Table of monocular right eye saccades

% Handle cases where one or both eyes have no saccades
if isempty(sacl) || isempty(sacr)
    msac = table();
    monol = sacl;
    monor = sacr;
    return;
end

% 1. Determine the temporal clusters (times where either eye is moving)
maxTime = max([max(sacl.offsetSample), max(sacr.offsetSample)]);
%make s, a vector of time with 1s whenever a microsaccade was happening in either eye
s = zeros(1, maxTime + 1);

for i = 1:height(sacl)
    s(sacl.onsetSample(i):sacl.offsetSample(i)) = 1;
end
for i = 1:height(sacr)
    s(sacr.onsetSample(i):sacr.offsetSample(i)) = 1;
end

% Find cluster onsets and offsets
m = find(diff([0, s, 0])); %Why add zeros to start and end?
m = reshape(m, 2, [])';     %reshape m so that onsets are in column 1 and offsets in column 2
m(:, 2) = m(:, 2) - 1; % Adjust for diff offset
numClusters = size(m, 1);

tempBinoc = [];
monolIdx = true(height(sacl), 1);
monorIdx = true(height(sacr), 1);

% 2. Identify Binocular vs Monocular clusters
for i = 1:numClusters
    % Find which saccades from each eye fall within this cluster
    rRows = find(sacr.onsetSample >= m(i,1) & sacr.offsetSample <= m(i,2));
    lRows = find(sacl.onsetSample >= m(i,1) & sacl.offsetSample <= m(i,2));

    if ~isempty(rRows) && ~isempty(lRows)
        % It's binocular: find the "main" saccade in each eye (by amplitude)
        [~, ir] = max(sacr.amp(rRows));
        [~, il] = max(sacl.amp(lRows));

        % Create combined row: Right eye data followed by Left eye data
        % Use cluster bounds for onset/offset to ensure they match
        rowR = sacr(rRows(ir), :);
        rowL = sacl(lRows(il), :);

        % Update onset/offset to the cluster boundaries
        rowR.onsetSample = m(i,1); rowR.offsetSample = m(i,2);
        rowL.onsetSample = m(i,1); rowL.offsetSample = m(i,2);
        %r and l onsets and offsets are the same, and cover the whole
        %cluster, so onset is the earliest time either eye started
        %moving and offset is the latest time either eye stopped.



        %this row now has the combined data and the L and R data saved in
        %it too

        tempBinoc = [tempBinoc; renameVars(rowR, 'R_'), renameVars(rowL, 'L_')];

        % Mark these as 'not monocular'
        monorIdx(rRows) = false;
        monolIdx(lRows) = false;
    end
end

% Assign monocular outputs
monol = sacl(monolIdx, :);
monor = sacr(monorIdx, :);

% 3. Final Merge of binocular events based on mergeInt
if ~isempty(tempBinoc)
    msac = tempBinoc(1, :);
    curr = 1;
    for i = 2:height(tempBinoc)
        % Check if gap between binocular events is small enough
        if tempBinoc.R_onsetSample(i) - msac.R_offsetSample(curr) <= mergeInt
            % Merge: update offsets
            msac.R_offsetSample(curr) = tempBinoc.R_offsetSample(i);
            msac.L_offsetSample(curr) = tempBinoc.L_offsetSample(i);
            % (Optionally update peak velocity/amplitude here if needed)
        else
            curr = curr + 1;
            msac(curr, :) = tempBinoc(i, :);
        end
    end

    % Aveage over eyes, to have 1 set of parameters for each binocular sacacde
    varsToAvg = {'onsetSample','offsetSample','peakVelocity','startX','startY','endX','endY','dx','dy','amp','totalAmpX','totalAmpY','maxCurveDeviation','curveRatio'};
    for vi=length(varsToAvg):(-1):1
        var = varsToAvg{vi};
        eval(sprintf('msac.%s = mean([msac.L_%s msac.R_%s],2);', var, var, var));
        %move this one to the front of table
        try
            eval(sprintf('msac = movevars(msac, ''%s'', ''Before'', ''%s'');', var, msac.Properties.VariableNames{1}));
        catch
            keyboard
        end
    end

    %add duration of each sacade
    msac.dur = msac.offsetSample-msac.onsetSample+1;
    msac = movevars(msac, "dur", 'Before',"peakVelocity");

    %add angle of each saccade, based on averaged change in x and
    %y-position
    msac.angle = 180/pi*atan2(msac.dy,msac.dx);
    msac = movevars(msac, "angle", 'Before',"totalAmpX");

else
    msac = table();
end

end

function T = renameVars(T, prefix)
T.Properties.VariableNames = strcat(prefix, T.Properties.VariableNames);
end
function [iteration] = funcReadContinue()
% Reads the iteration value of the existing continue.xml file in the instance of restarting.

fid = fopen('continue.xml','r');
if fid == -1
    error('ERROR: cannot find continue.xml for restart. Did you intent to restart?');
end

while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'iteration');
    if lineIndex
        iteration = strsplit(line);
        iteration = str2double(iteration{3});
        break
    else
        % Not correct line, keep searching.
    end
end

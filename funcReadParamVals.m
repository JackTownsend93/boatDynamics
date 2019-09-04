function [characteristicLength, resolution, uLB, uRef, outIter, cpIter] = funcReadParamVals(paramsFilename)
% Search through .xml file to find important values for calculating iteration steps etc.

fid = fopen(paramsFilename,'r');
if fid == -1
    error('ERROR: cannot find params_template.xml.');
end

% Search for parameter name then extract the relevant number or string
% following that regexp. Note format should be:
%
%   <varName> varVal <varName>
%
% Note the spaces between the equals are important as the third space
% delimited item is taken.

% characteristicLength.
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'<characteristicLength>');
    if lineIndex
        characteristicLength = strsplit(line);
        characteristicLength = str2double(characteristicLength{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% resolution.
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'<resolution>');
    if lineIndex
        resolution = strsplit(line);
        resolution = str2double(resolution{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% uLB.
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'<uLB>');
    if lineIndex
        uLB = strsplit(line);
        uLB = str2double(uLB{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% uRef.
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'<uRef>');
    if lineIndex
        uRef = strsplit(line);
        uRef = str2double(uRef{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% outIter.
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'<outIter>');
    if lineIndex
        outIter = strsplit(line);
        outIter = str2double(outIter{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% cpIter.
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'<cpIter>');
    if lineIndex
        cpIter = strsplit(line);
        cpIter = str2double(cpIter{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);


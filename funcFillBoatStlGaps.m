function [] = funcFillBoatStlGaps()
% Boat geometry .stl files are only created and save when boat movement 
% occurs. This can leave gaps in the output files where pressures .vtk 
% files or interface .stl files do not have a corresponding boat.stl file.
% This function fixes that as a post-processing step.
%     1. Read from params.xml the frequency of Palabos output to disk and 
%        boat.stl file name.
%     2. Loop through all iteration numbers, noting where matching boat.stl
%        files occur and where spaces are.
%     3. For every iteration number a boat.stl file does not exist for, 
%        copy to a file of that iteration number name the previous existing
%        boat.stl.

% 1. Read from params.xml.
fid = fopen('params.xml');
if fid == -1
	error('ERROR(funcFillBoatStlGaps.m): params.xml not found.');
end

% Get boat.stl name and frequency of Palabos output to disk.
while ~feof(fid)
	line = fgetl(fid);
	lineBoatName = regexp(line,'boatStl');
	lineOutIter  = regexp(line,'outIter');
	if lineBoatName
		boatName = strsplit(line);
		boatName = boatName{3};
	elseif lineOutIter
		outIter = strsplit(line);
		outIter = str2double(outIter{3});
	else
		% If neither found, do nothing.
	end
end

% Remove .stl extension from boatName.
boatName = boatName(1:length(boatName)-4);

% Count number of output files in /tmp.
interfaceFiles = dir('tmp/interface_*');
numPalabosOutputs = (length({interfaceFiles.name}));

% List boat.stl files in /tmp.
boatStlFiles = dir(['tmp/',boatName,'_*']);
boatStlFileNames = struct2cell(boatStlFiles);
boatStlFileNames = boatStlFileNames(1,:)';

% 2. Loop through all iterations checking if a corresponding boat.stl file
% exists.
matchingBoatStl = zeros(numPalabosOutputs,1); % Pre-allocating.
for i = 1:numPalabosOutputs
    palabosIterNum = i*outIter-outIter;
    palabosIterNumPadded = num2str(palabosIterNum,'%08.0f');
    matchedIndex = regexp(boatStlFileNames,palabosIterNumPadded);
    if isempty(cell2mat(matchedIndex))
        % No boat.stl file exists for this palabos iteration - index.
        matchingBoatStl(i,1) = false;
    else
        % A boat.stl files already exists for this iteration - index.
        matchingBoatStl(i,1) = true;
    end
end

% 3. Fill in non-matching cases by copying forward the previous boat.stl file.
% Note: a special exception will be to be made for the first gap.
previousMatch = boatName; j = 1;
for i = 1:length(matchingBoatStl)
    palabosIterNum = i*outIter-outIter;
    palabosIterNumPadded = num2str(palabosIterNum,'%08.0f');
    if matchingBoatStl(i,1)
        % Match - set the corresponding filename as the new latest match.
        previousMatch = cell2mat(boatStlFileNames(j,1));
        j = j + 1;
    else
        % No match - copy previousMatch file to a new boat.stl name with
        % appropriate iteration number.
        if strcmp(previousMatch,boatName)
            % Special case, copy main boat.stl file to iternumbered file.
            system(['cp ',boatName,'.stl tmp/',boatName,'_',palabosIterNumPadded,'.stl']);
        else
            system(['cp tmp/',previousMatch,' tmp/',boatName,'_',palabosIterNumPadded,'.stl']);
        end
    end
end

function [restartIter,restartTime] = funcReadLastSim(restarting)
% Reads the iteration value of the existing continue.xml file in the instance of restarting.

% Find latest sim-xxxxxxxx.out file in /tmp and open for reading.
simFiles = dir('tmp/sim-*');
simFiles = struct2cell(simFiles);
% Check if sim files exist in case of restarting.
if restarting && isempty(simFiles)
	fprintf('%s: No previous sim data found, check "restarting" value.\n',mfilename);
end
% Pick out restart iteration from sim files.
numSimFiles = size(simFiles);
simFilesLast = simFiles{1,numSimFiles(2)};
restartIter = str2num(cell2mat(regexp(simFilesLast,'\d*','Match')));

% Open file.
fid = fopen(sprintf('tmp/%s',simFilesLast),'r');
if fid == -1
    error('ERROR: cannot find latest sim file in tmp/');
end

% Skip to end of file.
offset = -750;
numLines = 10;
fseek(fid,offset,'eof');
text = cell(numLines,1);
for i = 1:numLines
	text(i) = {fgetl(fid)};
end

% Concatenate into one string.
textCat = '';
for i = 1: numLines
    textCat = strcat(textCat,text(i));
end
textCat = textCat{1};

% Search through to find 't = ' pattern and identify the time 
% following it using regexp.
restartTime = str2double(regexp(textCat,'(?<= t = [^0-9]*)[0-9]*\.?[0-9]+', 'match'));
% regexp usage explained:
% (?<= t = [^0-9]*) : start matching after 't = ' followed by non-number
% characters.
% [0-9]*            : then match 0 or more characters in 0-9 range.
% \.?               : allow for an optional decimal point.
% [0-9]+            : catch the number characters after the decimal.

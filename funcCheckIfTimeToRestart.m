function tSim = funcCheckIfTimeToRestart(jobID, tNext)
% Function that determines what time the currently running job has
% reached. Returns t_current, which should be checked in the main 
% script (t_current > t_next --> time to restart).

% Get the name of the output file using the current job ID.
jobIDstring = num2str(jobID);
filename = strcat('slurm-',jobIDstring,'.out');

% Only need to read a chunk (the last ten lines or so) of the output 
% or will spend a long time reading big files. Open the output
% file and start reading offset by a number of bytes.
fid    = fopen(filename,'r');
offset = -500;          % Number of bytes by which to offset.
fseek(fid,offset,'eof');

% Create an empty cell and populate with the output chunk.
numLines = 10;
textSample = cell(numLines,1);
for i = 1:numLines
    textSample(i) = {fgetl(fid)}; 
end

% Concatenate sampled text into one string.
textSampleCat = '';
for i = 1: numLines
    textSampleCat = strcat(textSampleCat,textSample(i));
end
textSampleCat = textSampleCat{1};

% Search through to find 't = ' pattern and identify the time 
% following it using regexp.
tSim = str2double(regexp(textSampleCat,'(?<= t = [^0-9]*)[0-9]*\.?[0-9]+', 'match'));
% regexp usage explained:
% (?<= t = [^0-9]*) : start matching after 't = ' followed by non-number
% characters.
% [0-9]*            : then match 0 or more characters in 0-9 range.
% \.?               : allow for an optional decimal point.
% [0-9]+            : catch the number characters after the decimal.

% Print results and handle errors.
format longg
if isempty(tSim)
    fprintf('%s: Cannot read a time value from job %s.\n', mfilename, jobID);
    fprintf('%s: Job has either failed or is still initialising.', mfilename);
    fprintf('\n');
else
    fprintf('%s: Waiting to reach next restart time...\n', mfilename);
    fprintf('%s:         Current sim time = %fs.\n', mfilename, tSim);
    fprintf('%s:          Next checkpoint = %fs.\n', mfilename, tNext);
    fprintf('\n');
end

fclose(fid);

end
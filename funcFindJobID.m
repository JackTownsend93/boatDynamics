function jobID = funcFindJobID(cmdout)

% Finds the job ID of the job most recently submitted slurm job.
% Note: only works if slurm is the job scheduler of choice.

% Convert sacct output to a string and find the job ID number following the known pattern 'Submitted batch job'.
slurmOutputString = strtrim(cmdout);
searchPattern = '.*Submitted batch job ([0-9]+).*';
% matchedTokens should be a single entry cell array containing the job ID.
matchedTokens = regexp(slurmOutputString, searchPattern, 'tokens', 'once');

% Check that matchedTokens is valid (i.e.: not empty).
if isempty(matchedTokens)
    jobID = '';
    fprintf('%s: Could not find job ID from %s.\n', mfilename, slurmOutputString);
else
    jobID = matchedTokens{1};
    fprintf('%s: Job ID %s was found from %s.\n', mfilename, jobID, slurmOutputString);
end

end
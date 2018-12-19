function funcConfirmJobIsKilled(jobID)
% Function to determine if checkpointing is finished by checking if the current job is dead.

% There are three possibile returns to cmdout for the command: "squeue -j <jobID>":
% 1. Job is still running:
%      example: 
%      JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
%      37417   compute dynamics s.659413  R       6:38      2 scs[0047,0065]
%
% 2. Job has recently been killed:
%      example:
%      JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
%
% 3. Job has been killed for a while:
%      example:
%      slurm_load_jobs error: Invalid job id specified
%
% Determine if script is still running by searching for the regexp of the jobID.

% Request job information.
command = strcat('squeue -j ',num2str(jobID));
[status,cmdout] = system(command);

isKilled = false;

jobStillRunning = true;
while jobStillRunning
    jobStatus = regexp(cmdout,jobID);
    % If jobID can be found in queue, job is still running.
    if jobStatus
        fprintf('%s: Still checkpointing... \n', mfilename);
        pause(5);
        [status,cmdout] = system(command);
    % If not, job is dead.
    else
        jobStillRunning = false; % Redundant.
        fprintf('%s: Checkpointing finished.\n', mfilename);
            break
    end
end

end
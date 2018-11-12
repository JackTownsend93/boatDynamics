function isKilled = funcConfirmJobIsKilled(jobID)
% Function will confirm if the checkpointing process has finished
% and the job in the cluster queue has been killed.

% If isKilled is true, the checkpointing has finished.

% Request specific job information.
jobID = num2str(jobID);
% strcat ignores whitespace, so must be concatenated as cell.
command = strcat({'squeue -j'},{' '},{jobID});
[status,cmdout] = system(command{1});

isKilled = false;

while isKilled == false
    if length(cmdout) > 100
        fprintf('%s: waiting for checkpointing to finish...\n', mfilename);
        pause(5)
        [status,cmdout] = system(command{1});
    else
        isKilled == true;
        fprintf('%s: checkpointing finished.\n', mfilename);
        break
    end
end

% Now clean up the old checkpoint files. IN ANOTHER FUNCTION




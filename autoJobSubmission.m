%% Script to automate boatHullFormSolver submissions.

% Notes for user:
% - You need two batchfiles: 'batchfile_init' (used only for first sim), and 'batchfile' (used for restarts).
% - The simulation will run for as many iterations as is stated in the params.xml file. A 'finished' file will then
%   be generated, prompting the matlab script to abort the simulation and restart under new conditions.

%% DEFINITIONS AND USER INPUT - CHECK THESE BEFORE EACH RUN!
paramfile = 'params.xml';
batchfile = 'batchfile_init.job'; % initial batchfile does not call continue.xml.

maxTime     = input('Enter the maximum physical time to run for:\n');
dt          = input('How often should the sim restart (seconds)?:\n');
numOfSims   = maxTime/dt-1;
tNext       = dt;

% Dynamics Variables.
% Boat CofG:
CofG = [-6.0;   % x.
         0.8;   % y.
         0.0;]; % z.
% STL filename without extension (Note: will be overwritten when moved; save a backup STL file in the desired initial position).
filenameSTL = 'dv15';

%% USER OPTION - Variable speed?.
% changeSpeed = input('Would you like to vary inlet velocity across sims? (y/n)\n', 's');
% if changeSpeed == 'y'
%     fprintf('Number of sims about to be run = %f\n',numOfSims);
%     speedPrevious  = input('Enter the speed (m/s) in the initial params_init.xml file?\n'); %can be automated later.
%     speedIncrement = input('Enter the increment in speed between each restart:\n');
% else if changeSpeed =='n'
%     fprintf('No change in inlet velocity.\n');
% else 
%     error('ERROR: You must enter "y" or "n".');
%     end
% end

%% USER OPTION - Boat dynamics?
% boatDynamics = input('Enable boat dynamics? (y/n)\n', 's');
% if boatDynamics == 'y'
%     fprintf('Boat dynamics ENABLED.\n');
% else if boatDynamics == 'n';
%     fprintf('Boat dynamics DISABLED.\n');
% else
%     error('ERROR: You must enter "y" or "n".');
%     end
% end


%% SUBMIT INITIAL SIMULATION.
[status,cmdout] = system('sbatch batchfile_init.job');
% Find job ID from cmdout.
jobID = funcFindJobID(cmdout);
slurmOutName = strcat('slurm-',jobID,'.out');

%batchNum = 1; % Part of the temporary batchfile fudge.

%% SUBMISSION LOOP.
% Begin loop based on counting the number of sim restarts.
i = 0; % Initiate at zero restarts.
while i <= numOfSims 
    % Check that the initial sim is running and not queued.
    if exist(slurmOutName,'file') == 2
        % Periodically check the time step to see when it exceeds tNext.
        tSim = funcCheckIfTimeToRestart(jobID, tNext);
        if tSim >= tNext
            % Begin restart: start checkpointing process.
            system('rm checkpoint_*'); % Delete old checkpoint files first.
            system('touch abort');     % Trigger checkpointing by creating abort file.
            
            % Wait for checkpoint files to finish writing.
            funcConfirmJobIsKilled(jobID);
            
            % INSERT DYNAMICS AND STL MANIPULATION HERE.
            % 1. Dynamics determines the rotation and z-displacement.
	    %    [SET MANUALLY, FOR NOW]:
	         rotation_deg = 0.50*cosd(10*i);
            zDisplacement_m   = 0.08*sind(10*i);
	    % Tranform current STL CofG to origin.
            CofG = funcManipulateSTL(filenameSTL, CofG, rotation_deg, zDisplacement_m, i);
	    
	    % Move the old slurm.out file to /tmp to avoid cluttering the working directory.
	    system(['mv slurm-',jobID,'.out tmp/']);
            % Submit the continue batchfile and get new jobID.
            [status,cmdout] = system('sbatch batchfile.job');
            jobID = funcFindJobID(cmdout);
	    slurmOutName = strcat('slurm-',jobID,'.out');
            
            % BEGIN TEMPORARY BATCHFILE FUDGE: batchfile for each restart, batchfile_1, batchfile_2... etc.
            %batchCommand = strcat('sbatch batchfile_',num2str(batchNum),'.job');
            %[status,cmdout] = system(batchCommand);
            %jobID = funcFindJobID(cmdout);
            %slurmOutName = strcat('slurm-',jobID,'.out');
            %batchNum = batchNum + 1;
            % END TEMPORARY BATCHFILE FUDGE.

            tNext = tNext + dt; % Increment restart time.
            i = i + 1;          % Increment number of restarts made.
        else
            pause(5)
        end    
    else
        fprintf('Waiting for output file to be written. Job may still be queuing...\n');
        pause(30)
    end
end

%% FINALISE.
% Kill final job.
system(['scancel ',jobID]);
% Move the final slurm.out file to /tmp.
system(['mv slurm-',jobID,'.out tmp/']);

fprintf('\n\n');
fprintf('           %%-------------------------------------------%%\n');
fprintf('           The auto submission process is COMPLETE:-\n');
fprintf('                        Total time: %f seconds.\n', maxTime);
fprintf('             Time between restarts: %f seconds.\n', dt);
fprintf('              Total number of sims: %d.\n', numOfSims);
fprintf('                      Final job ID: %s.\n', jobID);
fprintf('           %%-------------------------------------------%%\n');
fprintf('\n\n');

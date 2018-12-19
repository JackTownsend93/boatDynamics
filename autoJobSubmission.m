%% Script to automate boatHullFormSolver submissions.

% Notes for user:
% - You need two batchfiles: 'batchfile_initi.job' (only for first sim), and 'batchfile.job' (for restarts).
% - The hydro simulation will run for as many iterations as is stated in the params.xml file. A 'finished' file will then be generated, prompting the matlab script to abort the simulation and restart under new conditions.

%% DEFINITIONS AND USER INPUT - CHECK THESE BEFORE EACH RUN!
  paramfile = 'params.xml';         % Input parameters file.
  batchfile = 'batchfile_init.job'; % Initial batchfile (which does not call "continue.xml" as an input arg).
    maxTime = 10.0;    % Length of overall dynamic sim.
         dt =  0.05;   % Time of each component hydro sim.
filenameSTL = 'dv15';  % Also save a filename_init.stl file, as this one will be overwritten.
       CofG = [-6.0;   % x.
                0.8;   % y.
                0.0;]; % z.

% Calculate number of submits.
numOfSims = maxTime/dt-1;
    tNext = dt;

%% USER OPTION: Boat dynamics?
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

%% SUBMISSION LOOP.
% Begin loop based on counting the number of sim restarts.
i = 0; % Initiate at zero restarts.
while i <= numOfSims 
    % Check that the initial sim is running and not queued.
    if exist(slurmOutName,'file') == 2
    fprintf('Waiting for output file to be written. Job may still be queuing...\n');
        % Periodically check the time step to see when it exceeds tNext.
        tSim = funcCheckIfTimeToRestart(jobID, tNext);
        if tSim >= tNext
            % Begin restart: start checkpointing process.
            system('rm checkpoint_*'); % Delete old checkpoint files first.
            system('touch abort');     % Trigger checkpointing by creating abort file.
            
            % Wait for checkpoint files to finish writing.
            funcConfirmJobIsKilled(jobID);
            
            % DYNAMICS.
            % Dynamics function will determine the rotation and z-displacement.
	    % [SET MANUALLY, FOR NOW]:
	         rotation_deg = 0.50*cosd(10*i);
            zDisplacement_m   = 0.08*sind(10*i);
	    % Tranform current STL CofG to origin.
            CofG = funcManipulateSTL(filenameSTL, CofG, rotation_deg, zDisplacement_m);
	    
	    % Move the old slurm .out and .err files to /tmp to avoid cluttering the working directory.
	    system(['mv slurm-',jobID,'.out tmp/']);
	    system(['mv slurm-',jobID,'.err tmp/']);
            % Submit the continue batchfile and get new jobID.
            [status,cmdout] = system('sbatch batchfile.job');
            jobID = funcFindJobID(cmdout);
	    slurmOutName = strcat('slurm-',jobID,'.out');
            
            tNext = tNext + dt; % Increment restart time.
            i = i + 1;          % Increment number of restarts made.
        else
            pause(5)
        end    
    else
        pause(30)
    end
end

%% FINALISE.
% Kill final job.
system(['scancel ',jobID]);
% Move the final slurm.out file to /tmp.
system(['mv slurm-',jobID,'.out tmp/']);
% Replace moved .stl file with initial .stl file for next run.
system(['cp ',filenameSTL,'_init.stl ',filename,'.stl']);

fprintf('\n\n');
fprintf('           %%-------------------------------------------%%\n');
fprintf('           The auto submission process is COMPLETE:-\n');
fprintf('                        Total time: %f seconds.\n', maxTime);
fprintf('             Time between restarts: %f seconds.\n', dt);
fprintf('              Total number of sims: %d.\n', numOfSims);
fprintf('                      Final job ID: %s.\n', jobID);
fprintf('           %%-------------------------------------------%%\n');
fprintf('\n\n');

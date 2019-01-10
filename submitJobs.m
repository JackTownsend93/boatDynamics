%% Script to automate boatHullFormSolver submissions.

% ARE YOU RESTARTING?
restarting = false;

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

max_iter = 100;

%% SUBMISSION LOOP.
for i = 1: numOfSims 
	% Remove old params.xml file and write new one using template.
	system('rm -f params.xml');
	system(''); % Write max_iter value to MAX_ITER of params_template and save as params.xml.

	% Execute solver.
	system('mpiexec ./boatHullFormSolver params.xml continue.xml'); % BUT DO NOT SUBMIT WITH continue.xml IF IT IS FIRST JOB OF A NON-RESTART SIM.
	
	% Begin restart: start checkpointing process.
	%system('rm checkpoint_*'); % Delete old checkpoint files first.
	%system('touch abort');     % Trigger checkpointing by creating abort file.

	% DYNAMICS.
	% Dynamics function will determine the rotation and z-displacement.
	% [SET MANUALLY, FOR NOW]:
	% rotation_deg = 0.50*cosd(10*i);
	% zDisplacement_m   = 0.08*sind(10*i);
	% Tranform current STL CofG to origin.
	% CofG = funcManipulateSTL(filenameSTL, CofG, rotation_deg, zDisplacement_m);
end

%% FINALISE.
% Replace moved .stl file with initial .stl file for next run.
system(['cp ',filenameSTL,'_init.stl ',filenameSTL,'.stl']);
% Fill in Palabos output iterations that do not have a corresponding boat.stl file.
% funcFillBoatStlGaps(); 

fprintf('\n\n');
fprintf('           %%-------------------------------------------%%\n');
fprintf('           The auto submission process is COMPLETE:-\n');
fprintf('                        Total time: %f seconds.\n', maxTime);
fprintf('             Time between restarts: %f seconds.\n', dt);
fprintf('              Total number of sims: %d.\n', numOfSims);
fprintf('           %%-------------------------------------------%%\n');
fprintf('\n\n');

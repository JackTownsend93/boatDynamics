%% Script to automate boatHullFormSolver submissions.

%% DEFINITIONS AND USER INPUT - CHECK THESE BEFORE EACH RUN!
% General.
filenameSTL = 'dv15';  % .stl filename, no file extension.
restarting = false;    % Restarting from a continue.xml file? NOTE: Use same maxTime and dt if so.

% Dynamics.
maxTime =  30.0;   % Length of overall dynamic sim.
     dt =  0.5;   % Time of each component hydro sim.
   CofG = [-6.0;   % x.
            0.8;   % y.
            0.0;]; % z.

%% CALCULATE FROM INPUTS.
% Calculate number of submits.
numOfSims = maxTime/dt;
% Read params_template for necessary values.
[charL, resolution, uLB, uRef, outIter, cpIter] = funcReadParamVals('params_template.xml');
% Calculate Palabos timesteps.
dx_palabos = charL/(resolution-1);
dt_palabos = uLB*dx_palabos/uRef;
% Calculate number of iterations required per Palabos sim to reach dynamics dt for each sim.
max_iter_per_sim = round(dt/dt_palabos); % Each Palabos sim will run for this many iters to reach user-set dynamics dt.
stat_iter = round(max_iter_per_sim/100); % Rate of Palabos iteration status write to output file.
out_iter  = round(max_iter_per_sim/5);   % Rate of Palabos iteration full output to disk.
cp_iter   = max_iter_per_sim;            % Rate of checkpointing (once at the end of each sim).
% If restarting, fine the time at which the last sim finished.
if restarting
	[restartIter,restartTime] = funcReadLastSim();
else
	restartIter = 0; restartTime = 0;
end

j = 0; % Temporary preset motion variable.

%% SUBMISSION LOOP.
for i = 1+restartIter: numOfSims+restartIter
	% Calculate sim-specific iteration to terminate.
	max_iter_this_sim = max_iter_per_sim*i;

	% Stream values to params.xml file from params_template.xml file.
	system('rm -f params.xml');
	system(sprintf('sed ''s/MAX_ITER/%d/; s/STAT_ITER/%d/; s/OUT_ITER/%d/; s/CP_ITER/%d/'' params_template.xml > params.xml',max_iter_this_sim+1,stat_iter,out_iter,cp_iter)); % +1 to max_iter_this_sim to ensure checkpointing occurs correctly.
	
	system(sprintf('echo sim-%08d: running...',i));
	system(sprintf('echo sim-%08d: iterations: %d - %d.\n',i,max_iter_this_sim-max_iter_per_sim,max_iter_this_sim));
	system(sprintf('echo sim-%08d:       time: %.3f - %.3f s.\n',i,dt*i-dt,dt*i));

	% Execute solver.
	% On initial sim, do not use continue.xml as an input unless restarting from existing continue.xml file.
	if i == 1+restartIter && ~restarting;
		system(sprintf('mpiexec ./boatHullFormSolver params.xml > sim-%08d.out',i));
	else	
		system(sprintf('mpiexec ./boatHullFormSolver params.xml continue.xml > sim-%08d.out',i));
	end

	system(sprintf('echo sim-%08d: COMPLETE.',i));
	
	% Delete old checkpoint files (requires knowing last iter num).
	if i+restartIter > 1 % No checkpoint files on initial loop.
		% Get and format iteration number of sim i-1.
		iter_num_old = max_iter_this_sim-max_iter_per_sim;
		iter_num_old = num2str(iter_num_old,'%08.0f');
		% Use this to remove the old checkpoint file.
		system(sprintf('rm checkpoint_%s_*',iter_num_old));
	end
	
	% Move sim ouputs to /tmp.
	system(sprintf('mv sim-%08d.out tmp/',i));
	
	% DYNAMICS.
	% Dynamics function will determine the rotation and z-displacement.
	% [SET MANUALLY, FOR NOW]:
	j = j+1;
	rotation_deg = 0.50*cosd(10*j);
	zDisplacement_m = 0.05*sind(10*j);
	% Tranform current STL CofG to origin.
	CofG = funcManipulateSTL(filenameSTL, CofG, rotation_deg, zDisplacement_m);
end

%% CLEAN UP.
% Replace moved .stl file with initial .stl file for next run.
system(['cp ',filenameSTL,'_init.stl ',filenameSTL,'.stl']);
% Fill in Palabos output iterations that do not have a corresponding boat.stl file.
funcFillBoatStlGaps(); 
% Move the params.xml file to /tmp.
system('mv params.xml tmp/');

fprintf('\n\n');
fprintf('           %%-------------------------------------------%%\n');
fprintf('           The auto submission process is COMPLETE:-\n');
fprintf('           \n');
fprintf('                          Geometry: %s.\n', filenameSTL);
fprintf('           \n');
fprintf('                     Starting time: %f s.\n', restartTime);
fprintf('                    Finishing time: %f s.\n', restartTime+maxTime);
fprintf('                        Total time: %f s.\n', maxTime);
fprintf('             Time between restarts: %f s.\n', dt);
fprintf('              Total number of sims: %d.\n', numOfSims);
fprintf('           \n');
fprintf('           Palabos Stats:\n');
fprintf('                        Palabos dt: %f s.\n', dt_palabos);
fprintf('                        Palabos dx: %f s.\n', dx_palabos);
fprintf('                        Resolution: %d.\n', resolution);
fprintf('             Characteristic length: %f m.\n', charL);
fprintf('                     Lattice speed: %f m/s.\n', uLB);
fprintf('                   Reference speed: %f m/s.\n', uRef);
fprintf('            Rate of output to disk: %d iterations.\n', out_iter);
fprintf('           \n');
fprintf('           \n');
fprintf('           %%-------------------------------------------%%\n');
fprintf('\n\n');

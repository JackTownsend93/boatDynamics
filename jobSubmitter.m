% Rigid Body Dynamics for Palabos boatHullFormSolver.
% Jack Townsend - EngD Project.

% Rigid body dynamics submission script for use with Palabos boatHullFormSolver.exe. This script and 
% associated functions automates the submission of boatHullFormSolver jobs to the cluster wrapped in a 
% forward Euler dynamics solution.

% Instructions:
% The files the user will need to edit:
%     - <geometry_name>.stl
%         Have this in the working directory and ensure the .stl extension is in lower case. Also 
%         create a copy of this file of the format "<geometry_name>_init.stl". At the end of each 
%         dynamics sim the moved <geometry_name>.stl file will be replaced with the _init.stl version
%         so that the next sim will start with the boat in the intended initial position, not where the 
%         previous sim left off. Note that if the dynamics sim is killed before completion you will
%         need to do this manually (simply cp <geometry_name>_init.stl to <geometry_name>.stl. The 
%         _init.stl version is never actually used except to automate resetting of the boat to its 
%         original position.
%     - params_template.xml
%         This file configures all the paramters for the Palabos sims. Edit this file to your needs, 
%         leaving the "Iter" values as they are (they are calculated and streamed to a params.xml file
%         by this script.         
%     - batchfile.job
%         This is a standard batchfile for submitting this script to the cluster. One core is reserved 
%         for the running of this script, the solver is then executed on the remaining cores.
%     - /tmp directory
%         Ensure this directory exists before running or an error will be thrown. This contains the 
%         outputs of the simulation. Best practice is to clear this directory before submitting this 
%         script again. 

% Finally a few parameters in the next section of this script should be set by the user:
%   These are mostly self-explanatory. Note that maxTime and dt refer to the Euler dynamics solver, not 
%   the individual Palabos sims (those values are set in the params_template.xml file).
%   When the final Palabos sim has run, a number of checkpoint_xxxxxxxx_x.plb/.dat files and a 
%   continue.xml file will be left in the working directory. You can restart the dynamics sim from these 
%   by setting the restarting option below to true. Note that the maxTime and dt values must remain the 
%   as the previous sim.

%%--------------------- DEFINITIONS AND USER INPUT - CHECK THESE BEFORE EACH RUN! --------------------%%
% General.
filenameSTL = 'dv15';  % .stl filename, no file extension.
restarting = true;    % Restarting from a continue.xml file? (NOTE: Use same maxTime and dt if so).

% Dynamics parameters.
maxTime =  5.0;    % Length of overall dynamic sim.
     dt =  0.005;  % Time of each component hydro sim.
   CofG = [-5.0;   % x.
            0.8;   % y.
            0.0;]; % z.
m_boat = 10000;    % Mass of boat (kg).
   Iyy = 1;

%%--------------------------------------- CALCULATE FROM INPUTS --------------------------------------%%
% Calculate number of submits.
numOfSims = maxTime/dt;

% Read params_template for necessary values.
[charL, resolution, uLB, uRef, outIter, cpIter] = funcReadParamVals('params_template.xml');

% Calculate Palabos timesteps.
dx_palabos = charL/(resolution-1);
dt_palabos = uLB*dx_palabos/uRef;

% Calculate number of iterations required per Palabos sim to reach dynamics dt for each sim.
max_iter_per_sim = round(dt/dt_palabos); % Each Palabos sim will run for this many iters to reach user-set dynamics dt.
stat_iter = round(max_iter_per_sim/100)+1; % Rate of Palabos iteration status write to output file.
out_iter  = round(max_iter_per_sim/5);   % Rate of Palabos iteration full output to disk.
cp_iter   = max_iter_per_sim;            % Rate of checkpointing (once at the end of each sim).

% If restarting, find the time at which the last sim finished.
if restarting
	[restartIter,restartTime] = funcReadLastSim(restarting);
else
	restartIter = 0; restartTime = 0;
end

% Initial values for dynamics.
    z = [0,0];     z_dot = [0,0];
theta = [0,0]; theta_dot = [0,0];

%%========================================= SUBMISSION LOOP ==========================================%%
for i = 1+restartIter : numOfSims+restartIter

	% Calculate sim-specific iteration to terminate.
	max_iter_this_sim = max_iter_per_sim*i;
        cp_iter = max_iter_this_sim;
        
	% Stream values to params.xml file from params_template.xml file.
	system('rm -f params.xml');
	system(sprintf('sed ''s/MAX_ITER/%d/; s/STAT_ITER/%d/; s/OUT_ITER/%d/; s/CP_ITER/%d/; s/YVEL/%d/;'' params_template.xml > params.xml',max_iter_this_sim+1,stat_iter,out_iter,cp_iter,z_dot(2) )); % +1 to max_iter_this_sim to ensure checkpointing occurs correctly.

        system(sprintf('echo --------------------------------------------------------------------------------\n'));
	system(sprintf('echo ITERATION %d:\n',i));
        system(sprintf('echo sim-%08d: running...',i));
	system(sprintf('echo sim-%08d: iterations: %d - %d\n',i,max_iter_this_sim-max_iter_per_sim,max_iter_this_sim));
	system(sprintf('echo sim-%08d:       time: %.3f - %.3f s\n',i,dt*i-dt,dt*i));

	% Execute solver.
	% On initial sim, do not use continue.xml as an input unless restarting from existing continue.xml file.
	if i == 1+restartIter && ~restarting;
		system(sprintf('srun ./boatHullFormSolver params.xml > sim-%08d.out',i));
	else	
		system(sprintf('srun ./boatHullFormSolver params.xml continue.xml > sim-%08d.out',i));
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
	
	%%----------------------------------------- DYNAMICS -----------------------------------------%%
	% Dynamics function will determine the rotation and z-displacement.
	
        % %PRE-DETERMINED MOTION FOR TESTING.
	% j = j+1;
        % t = dt*(i-1);
        % rotation_deg    = 0.12*cosd(200*t);
        % zDisplacement_m = 0.03*sind(200*t);
        	
	% Read .stl vertices and vertex normals.
	[vertices, vertexNormals] = funcReadVertices(filenameSTL);
	% Read forces and pressures.
	[x_forceAvg, y_forceAvg, pressureCoords, pressureConnecs, pressureData] = funcReadForcesAndPressures(vertices,vertexNormals,max_iter_this_sim,out_iter);
	% Use dual mesh to acquire areas associated with pressures.
	% [pressureAreas, pressureVertexNormals] = funcDualMesh(pressureCoords, pressureConnecs, pressureData, filenameSTL);
	% Calculate moment.
	% [M_z] = funcCalculateMoment(CofG, pressureVertexNormals, pressureCoords, pressureData, pressureAreas);
        M_z = 0;

        % Determine boat motion.
        [theta_diff, theta, theta_dot, z_diff, z, z_dot] = funcMoveBoat(dt, y_forceAvg, M_z, m_boat, Iyy, z, z_dot, theta, theta_dot);
         
	% Tranform current STL CofG to origin.
	[CofG, iterationNum] = funcManipulateSTL(filenameSTL, CofG, theta_diff, z_diff);
        %----------------------------------------------------------------------------------------------%
end

%% ------------------------------------------ POST-PROCESSING ----------------------------------------%%
% Replace moved .stl file with initial .stl file for next run.
%system(['cp ',filenameSTL,'_init.stl ',filenameSTL,'.stl']);
% Fill in Palabos output iterations that do not have a corresponding boat.stl file.
funcFillBoatStlGaps(); 
% Move the params.xml file to /tmp.
system('mv params.xml tmp/');

fprintf('\n\n');
fprintf('           %%-------------------------------------------%%\n');
fprintf('           The auto submission process is COMPLETE:-\n');
fprintf('           \n');
fprintf('                          Geometry: %s\n', filenameSTL);
fprintf('           \n');
fprintf('                     Starting time: %f s\n', restartTime);
fprintf('                    Finishing time: %f s\n', restartTime+maxTime);
fprintf('                        Total time: %f s\n', maxTime);
fprintf('             Time between restarts: %f s\n', dt);
fprintf('              Total number of sims: %d\n', numOfSims);
fprintf('           \n');
fprintf('           Palabos Stats:\n');
fprintf('                        Palabos dt: %f s\n', dt_palabos);
fprintf('                        Palabos dx: %f s\n', dx_palabos);
fprintf('                        Resolution: %d\n', resolution);
fprintf('             Characteristic length: %f m\n', charL);
fprintf('                     Lattice speed: %f m/s\n', uLB);
fprintf('                   Reference speed: %f m/s\n', uRef);
fprintf('            Rate of output to disk: %d iterations\n', out_iter);
fprintf('           \n');
fprintf('           \n');
fprintf('           %%-------------------------------------------%%\n');
fprintf('\n\n');

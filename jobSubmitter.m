%% Script to automate boatHullFormSolver submissions.

% ARE YOU RESTARTING?
restarting = false;

%% DEFINITIONS AND USER INPUT - CHECK THESE BEFORE EACH RUN!
    maxTime = 10.0;    % Length of overall dynamic sim.
         dt =  0.05;   % Time of each component hydro sim.
filenameSTL = 'dv15';  % Also save a filename_init.stl file, as this one will be overwritten.
       CofG = [-6.0;   % x.
                0.8;   % y.
                0.0;]; % z.

% Calculate number of submits.
numOfSims = maxTime/dt-1;

% Read params_template for critical values.
[charL, resolution, uLB, uRef, outIter, cpIter] = funcReadParamVals('params_template.xml');

% Calculate number of iterations per Palabos simulation.
% dx_palabos = charL/(resolution-1);
% dt_palabos = uLB*dx_palabos/uRef;
max_iter = 500; % Calc this later, predefined for now.

%% SUBMISSION LOOP.
for i = 1: numOfSims
	%iterNum =  
	% Remove old params.xml file and write new one using template.
	system('rm -f params.xml');
	system(sprintf('sed ''s/MAX_ITER/%d/g'' params_template.xml > params.xml',max_iter));
	
	% Execute solver.
	% On initial sim, do not use continue.xml as an input arg.
	if i == 1 && ~restarting;
		system(sprintf('mpiexec ./boatHullFormSolver params.xml > sim-%d.out',i));
	else	
		system(sprintf('mpiexec ./boatHullFormSolver params.xml continue.xml > sim-%d.out',i));
	end
	
	% Delete old checkpoint files (requires knowing last iter num).
	%system('rm checkpoint_*'); % Delete old checkpoint files first.

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
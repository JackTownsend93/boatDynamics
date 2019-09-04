function [x_forceAvg, y_forceAvg, pressureCoords, pressureConnecs, pressureData] = funcReadForcesAndPressures(vertices,vertexNormals,max_iter_this_sim,out_iter)
%% Reads forces and pressures from Palabos outputs.
% Forces are contained in total_force_on_boat.dat file in /tmp.
% Pressures are contained in boat_pressure_xxxxxxxx.vtk (xxxxxxxx =
% iteration number of last output).

%% READ FORCES. 
fid = fopen('tmp/total_force_on_boat.dat');
if fid == -1
    error('%s: ERROR tmp/total_forces_on_boat.dat not found.',mfilename);
end
fprintf('%s: Reading forces...\n',mfilename);

% Determine how many outputs are made during this Palabos simulation.
num_outputs_this_sim = floor(max_iter_this_sim/out_iter);
% Take the last 10% of force values to be averaged.
avgRange = floor(num_outputs_this_sim*0.1);
% Rounded down, but must avoid zero val.
if avgRange == 0
    avgRange = 1;
end

% Count total lines in force file.
totalLines = 0;
while ~feof(fid)
    line = fgetl(fid);
    totalLines = totalLines + 1;
end
beginRange = totalLines - avgRange;

% Rewind and count through to averaging range, pick values to be averaged.
frewind(fid);
lineCount = 0; avgIndex = 1;
x_forcesToAvg = zeros(avgRange,1);
y_forcesToAvg = zeros(avgRange,1);
while ~feof(fid)
    if lineCount < beginRange
        line = fgetl(fid);
        lineCount = lineCount + 1;
    else
        line = strsplit(fgetl(fid));
        x_forcesToAvg(avgIndex,1) = str2double(line{1,3});
        y_forcesToAvg(avgIndex,1) = str2double(line{1,4});
        lineCount = lineCount + 1;
        avgIndex = avgIndex + 1;
    end
end

% Average the values.
x_forceAvg = sum(x_forcesToAvg)/length(x_forcesToAvg);
y_forceAvg = sum(y_forcesToAvg)/length(y_forcesToAvg);
fprintf('%s: The average DRAG for the last %d Palabos iterations is %fN.\n',mfilename,length(x_forcesToAvg),x_forceAvg);
fprintf('%s: The average LIFT for the last %d Palabos iterations is %fN.\n',mfilename,length(y_forcesToAvg),y_forceAvg);

%% READ PRESSURES.
% First need to identify which is the latest boat pressure .vtk file.
% Format is:    boat_pressure_xxxxxxxx.vtk
% Where xxxxxxxx is an eight-digit number of the last Palabos iteration at 
% which output files were written.

% Find latest pressure.vtk file.
pressureFiles = dir('tmp/boat_pressure_*');
filenamePressure = pressureFiles(length({pressureFiles.name})).name;

% Open this file to read pressure data.
fid = fopen(['tmp/',filenamePressure]);
if fid == -1
    error('%s: ERROR tmp/boat_pressure_xxxxxxxx.vtk not found.',mfilename);
end

% Read number of points.
regexpFound = [];
while isempty(regexpFound)
    line = fgetl(fid);
    regexpFound = regexp(line,'POINTS');
    if regexpFound
        % Found line for numPoints.
        numPoints = strsplit(line);
        numPoints = str2double(numPoints{2});
    end
end

% Read all point coords.
fprintf('%s: Reading pressure coords...\n',mfilename);
pressureCoords = zeros(numPoints,3);
for i = 1:numPoints
    line = strsplit(fgetl(fid));
    pressureCoords(i,:) = str2double(line);
end

% Skip a line.
fgetl(fid);

% Read number of connectivities.
line = fgetl(fid);
numConnecs = strsplit(line);
numConnecs = str2double(numConnecs{2});

% Read connectivities.
fprintf('%s: Reading connectivities...\n',mfilename);
pressureConnecs = zeros(numConnecs,4);
for i = 1:numConnecs
	line = strsplit(fgetl(fid));
	pressureConnecs(i,:) = str2double(line);
end
pressureConnecs = pressureConnecs(:,2:4); % First column is redundant.

% Read up to "POINT_DATA" and check that numPoint == numPressures.
regexpFound = [];
while isempty(regexpFound)
    line = fgetl(fid);
    regexpFound = regexp(line, 'POINT_DATA');
    if regexpFound
        % Found line for numPressures.
        numPressures = strsplit(line);
        numPressures = str2double(numPressures{2});
    end
end
if numPoints ~= numPressures
    error('%s: Number of points does not match number of pressures in boat_pressure.vtk file',mfilename);
end
% Skip two lines.
fgetl(fid); fgetl(fid);

fprintf('%s: Reading pressure data...\n',mfilename);
% Read all pressure data.
pressureData = zeros(numPoints,1);
for i = 1:numPoints
    line = strsplit(fgetl(fid));
    pressureData(i,:) = str2double(line);
    if mod(i,1000000) == 0
      fprintf('%s: 1,000,000 of %d pressures read.\n',mfilename,numPoints);
    end
end

% DEBUGGING: Plotting pressure at coords.
%scatter3(pressureCoords(:,1),pressureCoords(:,2),pressureCoords(:,3),1,pressureData(:,1));
%quiver3(vertices(:,1),vertices(:,2),vertices(:,3),vertexNormals(:,1),vertexNormals(:,2),vertexNormals(:,3));
%axis equal;
%hold on;

fclose(fid);

end

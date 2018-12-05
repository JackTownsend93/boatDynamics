function [x_forceAvg, y_forceAvg, pressureCoords, pressureData] = funcReadForcesAndPressures()
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

% Format is simply a Nx4 matrix of columns: time, x_force, y_force, z_force.
% Take average of the last 10% of the lines in the file for y_force.

% Count number of lines in force file.
totalLines = 0;
while ~feof(fid)
    line = fgetl(fid);
    totalLines = totalLines + 1;
end

% Take last 10% of lines to get an average of y_force.
lastLinesBegin = round(totalLines - 0.1*(totalLines));
frewind(fid);
countLines = 0; numLinesToAvg = 0;
x_forcesToAvg = zeros(totalLines-lastLinesBegin-1,1);
y_forcesToAvg = zeros(totalLines-lastLinesBegin-1,1);
while ~feof(fid)
    if countLines < lastLinesBegin
        line = fgetl(fid);
        countLines = countLines + 1;
    else 
        numLinesToAvg = numLinesToAvg + 1;
	line = strsplit(fgetl(fid));
        x_forcesToAvg(numLinesToAvg,1) = str2double(line{1,3});
        y_forcesToAvg(numLinesToAvg,1) = str2double(line{1,4});
    end
end

% Averaging.
x_forceAvg = sum(x_forcesToAvg)/length(x_forcesToAvg);
y_forceAvg = sum(y_forcesToAvg)/length(y_forcesToAvg);
fprintf('%s: The average DRAG for the last %d Palabos iterations is %fN.\n',mfilename,length(x_forcesToAvg),x_forceAvg);
fprintf('%s: The average LIFT for the last %d Palabos iterations is %fN.\n',mfilename,length(y_forcesToAvg),y_forceAvg);
fclose(fid);

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
end

% Plotting pressure at coords.
scatter3(pressureCoords(:,1),pressureCoords(:,2),pressureCoords(:,3),1,pressureData(:,1));
axis equal;

fclose(fid);



% UNUSED HERE NOW, BUT USEFUL; APPLY THIS TO POINTS WHERE FIXED PARAMS.XML LINES ARE
% ASSUMED.
% % First read user-set value to frequency of outputs.
% fid = fopen('params.xml');
% if fid == -1
%     error('%s: ERROR params.xml not found.',mfilename);
% end
% % Find outIter line and read value.
% while ~feof(fid)
%     line = fgetl(fid);
%     startIndex = regexp(line,'outIter');
%     if startIndex
%         % outIter line found, get value.
%         outIter = strsplit(line);
%         outIter = str2double(outIter{3});
%         break
%     else
%         % Not outIter line, continue searching.
%     end
% end

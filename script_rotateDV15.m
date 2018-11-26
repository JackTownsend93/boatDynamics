% Script to rotate an stl file.
clear; clc; close all;

% Define stl filename.
stlFilename = 'dv15.stl';

% Read in dv15.stl
[vertices, faces, normals, stlName] = stlRead(stlFilename);

% Define CofG.
CofG = [-6.0;   % x.
         0.8;   % y.
         0.0;]; % z.

% Manipulate geometry.
% User values:
        rotation_deg = 10.0;
vertDisplacement_m   =  0.0;

% Transform CofG to origin.
vertices(:,1) = vertices(:,1) - CofG(1);
vertices(:,2) = vertices(:,2) - CofG(2);
vertices(:,3) = vertices(:,3) - CofG(3);

% Rotation matrix for 3D rotation about z-axis.
Rz = [ cosd(rotation_deg), sind(rotation_deg), 0  ;
      -sind(rotation_deg), cosd(rotation_deg), 0  ;
                        0,                  0, 1 ];
% Apply rotation.
vertices = vertices * Rz';

% Revert transformation to pre-rotation position.
vertices(:,1) = vertices(:,1) + CofG(1);
vertices(:,2) = vertices(:,2) + CofG(2);
vertices(:,3) = vertices(:,3) + CofG(3);

% Apply vertical displacement.
vertices(:,2) = vertices(:,2)+vertDisplacement_m;

% Write to a new stl.
stlWrite('new.stl', faces, vertices);

% Plot new stl.
[verticesNew, facesNew, normalsNew, stlNameNew] = stlRead('new.stl');
stlPlot(verticesNew, facesNew, stlNameNew);
title(['Transformed by ', num2str(vertDisplacement_m) 'm in z-direction and rotated by ', num2str(rotation_deg), char(176), ' about z-axis.']);
xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis');
view([0  90]);

% Plot CofG.
hold on;
scatter3(CofG(1),CofG(2),CofG(3)+2,100,'rx')
legend('Boat Geometry','CofG');
set(gca,'XLim',[-15 3],'YLim',[-5 5],'ZLim',[-3 3])

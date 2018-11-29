function CofG = funcManipulateSTL(filenameSTL, CofG, rotation_deg, zDisplacement_m, i)
%% Manipulates the .stl geometry by rotating and transforming in z-direction according to dynamics.
%  1. Read in .stl file (binary or ASCII).
%  2. Move whole .stl so that the CofG is at the origin.
%  3. Perform rotation.
%  4. Tranform back to pre-rotation position.
%  5. Perform z translation.
%  6. Tranform CofG also.
%  7. Write to new .stl file and store a copy for post-processing later.

%% 1. Read .stl file.
[vertices, faces, normals, nameSTL] = stlRead([filenameSTL,'.stl']);

%% 2. Move whole .stl file so that the CofG is at the origin.
vertices(:,1) = vertices(:,1) - CofG(1);
vertices(:,2) = vertices(:,2) - CofG(2);
vertices(:,3) = vertices(:,3) - CofG(3);

%% 3. Perform rotation.
% Define rotation matrix for 3D rotation about the z-axis.
Rz = [ cosd(rotation_deg), sind(rotation_deg), 0  ;
      -sind(rotation_deg), cosd(rotation_deg), 0  ;
                        0,                  0, 1 ];
% Apply rotation to the vertices.
vertices = vertices * Rz';

%% 4. Return to pre-rotation position.
vertices(:,1) = vertices(:,1) + CofG(1);
vertices(:,2) = vertices(:,2) + CofG(2);
vertices(:,3) = vertices(:,3) + CofG(3);

%% 5. Z-translation.
vertices(:,2) = vertices(:,2)+zDisplacement_m;

%% 6. Transform CofG identically (no rotation required as rotation is effectively about CofG).
CofG(2) = CofG(2) + zDisplacement_m;

%% 7. Write .stl files.
% Write an updated .stl file out in working directory.
stlWrite([filenameSTL,'.stl'], faces, vertices);
% Store a copy in /stl_stored for paraview post-processing.
filenameSTL_i = [filenameSTL,'_',num2str(i)];
stlWrite(['tmp/',filenameSTL_i,'.stl'], faces, vertices);


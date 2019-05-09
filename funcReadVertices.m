function [vertices, vertexNormals] = funcReadVertices(filenameSTL)
%% Read vertices and their normals from the active .stl file.

% Read .stl file.
[vertices, faces, faceNormals, nameSTL] = stlRead([filenameSTL,'.stl']);

% Get vertex normals.
fprintf('%s: Reading .stl vertices and vertex normals.\n',mfilename);
vertexNormals = zeros(size(vertices));
% Add the normals for each point of every triangle.
vertexNormals(faces(:,1),:) = vertexNormals(faces(:,1),:) + faceNormals;
vertexNormals(faces(:,2),:) = vertexNormals(faces(:,2),:) + faceNormals;
vertexNormals(faces(:,3),:) = vertexNormals(faces(:,3),:) + faceNormals;
% Calculate vertex normals.
vertexNormals = vertexNormals./repmat(sqrt(sum(vertexNormals.^2,2)),[1,3]);
vertexNormals = -vertexNormals; % Invert to point inward of volume.

end

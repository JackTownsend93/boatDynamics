function [M] = funcCalculateMoment(CofG, vertexNormals, pressureCoords, pressureData)
% Calculates moment from pressures and .stl data.

mag = zeros(size(vertexNoramals));
mag(:,1) = vertexNormals(:,1)*pressureData';
mag(:,2) = vertexNormals(:,2)*pressureData';
mag(:,3) = vertexNormals(:,3)*pressureData';



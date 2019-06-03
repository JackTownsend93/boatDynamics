function [pressureAreas, pressureVertexNorms] = funcDualMesh(pressureCoords, pressureConnecs, pressureData)
%% Create a dual mesh to find areas associated with each pressure vector.

% Make the dual mesh.
pp = pressureCoords;
tt = pressureConnecs+1;
[cp,ce,dualNodalCoords,dualEdgeIndices] = makedual2(pp,tt);

% Calculate centroids and areas.
[dualCentroids,dualAreas] = geomdual2(cp,ce,dualNodalCoords,dualEdgeIndices);

% Find the nearest neighbour in dualCentorids for every query point in 
% pressureCoords, returning the indicesof the nearest neighbour 'idx'
% and the Euclidean distance 'd' (useful for verifications, ensure this 
% is small).
[idx,d] = knnsearch(dualCentroids,pressureCoords);
% Pre-allocate pressureAreas.
pressureAreas = zeros(length(pressureCoords),1);
% Cycle through idx and match.
for i = 1:length(idx)
    pressureAreas(i,1) = pressureAreas(i,1) + dualAreas(idx(i),1);
end
% Check conservation of area.
if abs(sum(pressureAreas) - sum(dualAreas)) < 0.0001
   fprintf('%s: Area correctly conserved.\n',mfilename);
else
   fprintf('%s: You have lost some dual cells somewhere!\n',mfilename);
end

end

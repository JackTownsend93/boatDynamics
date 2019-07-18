function [pressureAreas, pressureVertexNormals] = funcDualMesh(pressureCoords, pressureConnecs, pressureData, filenameSTL)
%% Create a dual mesh to find areas associated with each pressure vector.

tStart = tic; % Time the process.

%% Create the dual mesh.
% Rename for convention.
pp = pressureCoords;
tt = pressureConnecs +1; % Start from 1, not 0, due to Matlab indexing.

% Create the dual mesh.
[cp,ce,dualNodalCoords,dualEdgeIndices] = makedual2(pressureCoords,tt);
% Calculate areas and centroids of dual cells.
[dualCentroids,dualAreas] = geomdual2(cp,ce,dualNodalCoords,dualEdgeIndices);

%% Match .vtk pressures to dual mesh centroids.
% NOTE: # of dual mesh areas might != # of pressure coords points.
% KnnSearch avoiding duplicate associations (computationally slow).
% Make a copy of dualCentroids that will have rows eliminated.
tmp_dualCentroids = dualCentroids;
% Pre-alloc.
idx_dual = zeros(length(pressureCoords),1);
d_dual = zeros(length(pressureCoords),1);
for i = 1:length(pressureCoords)
    % Going through pressureCoords one by one, find the index of
    % the nearest dual centroid and store it in idx.
    [idx_dual(i,:),d_dual(i,:)] = knnsearch(tmp_dualCentroids,pressureCoords(i,:));
    % Now eliminate that tmp_dualCentroid row to avoid double 
    % associations.
    tmp_dualCentroids(idx_dual(i,:),:) = []; % Delete this row.
end
% Calc some factors to determine inaccuracy of associations.
% Average Euclidean distance between associations.
fprintf('Average distance between matches: %f \n',sum(d_dual)/length(d_dual));
% Maximum distance.
fprintf('Maximum distance between matches: %f \n',max(d_dual));

% Create a variable called pressureAreas, same length as pressure
% coords. Fill it using idx so that pressureAreas correspond to
% pressureCoords.
pressureAreas = zeros(length(pressureCoords),1);
for i = 1:length(idx_dual)
    pressureAreas(i,1) = pressureAreas(i,1) + dualAreas(idx_dual(i),1);
end

% Check conservation of area.
if abs(sum(pressureAreas) - sum(dualAreas)) < 0.0001
    fprintf('Area correctly conserved.\n');
else
    fprintf('You have lost some dual cells somewhere.\n');
    fprintf('                    Dual mesh total area: %f m^2\n',sum(dualAreas));
    fprintf('Total area associated to pressure values: %f m^2\n',sum(pressureAreas));
end

% For now, use the .stl normals instead of calculating .vtk normals. Use 
% knnsearch for this again, no need to avoid duplicating associations.
[vertices,vertexNormalsSTL] = funcReadVertices(filenameSTL);
[idx_norms,d_norms] = knnsearch(vertices,pressureCoords);
% Use idx_norms to match vertexNormals to correct position in
% pressureVertexNormals.
pressureVertexNormals = zeros(length(pressureCoords),3);
for i = 1:length(idx_norms)
    pressureVertexNormals(i,:) = -vertexNormalsSTL(idx_norms(i),:);
end

% Calc some factors to determine inaccuracy of associations.
% Average Euclidean distance between associations.
fprintf('Average distance between matches: %f \n',sum(d_norms)/length(d_norms));
% Maximum distance.
fprintf('Maximum distance between matches: %f \n',max(d_norms));

% Dual mesh timer.
tElapsed = toc(tStart);
tElapsedHours = tElapsed/3600;
tElapsedMinutes = tElapsed/60;
fprintf('%s: Time taken for dual mesh process: %dh:%dm:%ds.\n',mfilename,floor(tElapsedHours),floor(mod(tElapsedMinutes,60)),floor(mod(tElapsed,60)));

%% Plots for debugging (comment away when not debugging).
% Plot original triangulation.
%figure(1);
%hold on; title('Triangulation');
%drawtria2(pp,tt);
% % view(-150,20);
%set(gca,'units','normalized'); axis image off;

% Plot dual complex.
%figure(2);
%hold on; title('Dual complex');
%drawdual2(cp,ce,dualNodalCoords,dualEdgeIndices);
% % view(-150,20);
%set(gca,'units','normalized'); axis image off;

% Plot both overlaid.
%figure(3);
%hold on; title('Combine');
%drawtria2(pp,tt);
%set(gca,'units','normalized'); axis image off;
%hold on;
%drawdual2(cp,ce,dualNodalCoords,dualEdgeIndices);
%set(gca,'units','normalized'); axis image off;

% Plot cell centroids.
%figure(4);
%hold on; title('Cell centroids');
%drawdual2(cp,ce,dualNodalCoords,dualEdgeIndices);
%plot3(dualCentroids(:,1),dualCentroids(:,2),dualCentroids(:,3),'k.');
%set(gca,'units','normalized'); axis image off;

% Plot cell areas.
%figure(5); hold on;
%title('Cell areas');
%drawdual2(cp,ce,dualNodalCoords,dualEdgeIndices,dualAreas);
%set(gca,'units','normalized'); axis image off;

% Plot .vtk coordinates positions against dual area centroids.
%figure(6); hold on;
%title('Pressure .vtk Coords vs Dual Area Centroids');
%drawdual2(cp,ce,dualNodalCoords,dualEdgeIndices);
%plot3(dualCentroids(:,1),dualCentroids(:,2),dualCentroids(:,3),'k.');
%scatter3(pressureCoords(:,1),pressureCoords(:,2),pressureCoords(:,3),5,'rx');
%plot3(tmp_dualCentroids(:,1),tmp_dualCentroids(:,2),tmp_dualCentroids(:,3),'bx');
%legend('Dual Surface','Dual Cells','Dual Centroids','.vtk Pressure Coords','Unused centroids');
%set(gca,'units','normalized'); axis image off;

% Plot .stl vertex normals in position on the dual mesh.
%figure(7); hold on;
%title('.stl vertex normals in position on the dual mesh');
%drawdual2(cp,ce,dualNodalCoords,dualEdgeIndices);
%scatter3(pressureCoords(:,1),pressureCoords(:,2),pressureCoords(:,3),1,pressureData(:,1));
%quiver3(pressureCoords(:,1),pressureCoords(:,2),pressureCoords(:,3),pressureVertexNorms(:,1),pressureVertexNorms(:,2),pressureVertexNorms(:,3));
%set(gca,'units','normalized'); axis image off;

end

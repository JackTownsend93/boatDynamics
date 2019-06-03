function [M_z] = funcCalculateMoment(CofG, pressureVertexNormals, pressureCoords, pressureData, pressureAreas)
% Calculates the moment in about the z-axis.

% Calculate forces from pressure data.
f = pressureData.*pressureAreas;    % Scalar force at each point.
f_vec = [f,f,f].*pressureVertexNormals;     % Apply magnitude to unit normals.
f_res = f'*pressureVertexNormals;           % Calculate resulatant force components to validate against Palabos values.
fprintf('Drag:    %f N\nLift:    %f N\nSway:    %f N\n',f_res(1,1),f_res(1,2),f_res(1,3));
% quiver3(pressureCoords(:,1),pressureCoords(:,2),pressureCoords(:,3),f_vec(:,1),f_vec(:,2),f_vec(:,3)); % Visualise.

% Translate origin to CofG.
pressureCoords_CofG = pressureCoords; % Temporary.
pressureCoords_CofG(:,1) = pressureCoords(:,1)-CofG(1);
pressureCoords_CofG(:,2) = pressureCoords(:,2)-CofG(2);
pressureCoords_CofG(:,3) = pressureCoords(:,3)-CofG(3);
% quiver3(pressureCoords_CofG(:,1),pressureCoords_CofG(:,2),pressureCoords_CofG(:,3),f_vec(:,1),f_vec(:,2),f_vec(:,3)); % Visualise.

% Find pitching moment about CofG using moment arms given by pressureCoords_CofG.
M_z = sum(f_vec(:,1).*pressureCoords_CofG(:,2) + f_vec(:,2).*pressureCoords_CofG(:,1));

%% Plot verifying forces/moments.
%figure; hold on; grid on; axis equal;
%quiver3(pressureCoords_CofG(:,1),pressureCoords_CofG(:,2),pressureCoords_CofG(:,3),f_vec(:,1),f_vec(:,2),f_vec(:,3));
%scatter3(0,0,0,'rx');
%f_diag = diag(f_res);
%% Scale and plot.
%f_diag = f_diag./max(abs(f_res));
%quiver3([0,0,0],[0,0,0],[0,0,0],f_diag(1,:),f_diag(2,:),f_diag(3,:));

end

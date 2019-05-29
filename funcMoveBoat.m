function [rotation_deg, zDisplacement_m] = funcMoveBoat(dt, y_forceAvg, M_z, m_boat)
% Determine 2DoF boat motion (heave and pitch) given the global timestep and force/moment.

% Heave.
zDisplacement_m = (y_forceAvg-m_boat*9.81)/m_boat * dt^2;

% Pitch.
rotation_deg = 0; % Temporary.


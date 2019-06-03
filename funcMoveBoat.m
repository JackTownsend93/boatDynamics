function [theta_diff, theta, theta_dot, z_diff, z, z_dot] = funcMoveBoat(dt, y_forceAvg, M_z, m_boat, Iyy, z, z_dot, theta, theta_dot)
% Determine 2DoF boat motion (heave and pitch) given the global timestep and force/moment.

% Accel due to grav at sea level.
g = 9.81;

% Only need previous values of z and z_dot, so both these are 1x2 vectors where:
%     z(1) = z_i    and    z(2) = z_i+1
% And likewise for z_dot.

% Heave.
        F_z = y_forceAvg;
delta_z_dot = ((F_z-m_boat*g)/m_boat)*dt;
   z_dot(2) = z_dot(1) + delta_z_dot;
       z(2) = z(1) + z_dot(1)*dt;

% funcManipulateSTL.m operates about the boat's previous position, not global coords, so the 
% difference between z values is needed.
z_diff = z(2) - z(1);

% "Iterate" z values.
z_dot(1) = z_dot(2);
    z(1) = z(2);

%Pitch.
theta_diff = 0; % Temporary.
% 
% delta_theta_dot = (M_z/Iyy)*dt;
% theta_dot(2) = theta_dot(1) + delta_theta_dot;
% theta(2) = theta(1) + theta_dot(1)*dt;

% Boat is moved according to change in pitch, not absolute pitch.
% theta_diff = theta(2) - theta(1);

% "Iterating" theta values.
% theta_dot(1) = theta_dot(2);
%     theta(1) = theta(2);

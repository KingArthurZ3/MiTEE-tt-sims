clear; close all; clc;
global mu;
Re = 6371e3;            %earth radius
Ro = [Re+350e3 Re+250e3];          %orbit radius (350 km above sea level)
mu = 3.98576e14;        %Gravitational constant times earth mass (m^3 * s^-2)
w = [0 0;sqrt(mu/(Ro(1)^3)) sqrt(mu./(Ro(2)^3));0 0];    %initial orbital angular velocity
r_init = [0;0;1]*Ro;       %initial position in inertial frame (no y-component)
% currently defining inertial frame to be aligned with orbit (i.e.
% y-inertial is orthogonal to orbit plane). This is because we are not
% currently considering out-of-plane orbit (it makes the math really hard)
%make sure that whatever is in the vector before Ro is a unit vector!
% Ex: [0;0;-1].*Ro means the initial position is at the top of the paper


%rotate position 90 deg ccw about y-axis to get velocity direction (unit)
v_hat1 = [0 0 -1; 0 1 0; 1 0 0]*r_init(:,1) ./ norm([0 0 1; 0 1 0; -1 0 0]*r_init(:,1));
v_hat2 = [0 0 -1; 0 1 0; 1 0 0]*r_init(:,2) ./ norm([0 0 1; 0 1 0; -1 0 0]*r_init(:,2));
v_orbit1 = abs(w(2, 1).*Ro(1)); %initial linear orbital speed r1dot
v_orbit2 = abs(w(2, 2).*Ro(2)); %initial linear orbital speed r2dot
v_orbit_inertial1 = v_hat1 .* v_orbit1;    %initial linear orbital velocity vector (inertial frame)
v_orbit_inertial2 = v_hat2 .* v_orbit2;    %initial linear orbital velocity vector (inertial frame)

tspan = [0 5429];                       %orbital period is 5429 sec
X1_init = [r_init(:,1); v_orbit_inertial1];    % state vector
X2_init = [r_init(:,2); v_orbit_inertial2];    % state vector
[t_vec,X1] = ode45(@twobody, tspan, X1_init);
[t_vec,X2] = ode45(@twobody, tspan, X2_init);


figure(1);
hold on
scatter(X1(1:end,1),X1(1:end,3), 'r');
scatter(X2(1:end,1),X2(1:end,3), 'b');
axis equal;
pbaspect([1 1 1]);
hold off
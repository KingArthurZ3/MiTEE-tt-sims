
%Version 1.1
%Date modified: 2/6/2019 5:30 pm
%Notes:
%   -keep initial w small (~10^-2 rad/s)
%   -added nonhomogeneous solution, but for right now, set Mc = 0
%Tasks:
%   -Need to make simulation output attitude as well. Use the 1/cos matrix
%    in the AE347 notes
%   -Create model for moments and forces (e.g. control torque, gravity gradient,
%    Lorentz force etc.)
%   -implement torque controller
%   -Add Newton's second law to model
%   -Model and output state of each end body (i.e. attitude and angular
%    rate of picosat and nanosat)
%   -Calculate tension
%   -declare tether release controller and add to model
%   -implement tether release controller
%   -add documentation when possible
clear; close all; clc;
dt = 0.03;                  %simulation time  step (s)
method = 1;                 %1: Forward Euler, 2: RK4

% Spacecraft Specs to be implemented (main body is MiTEE 1 rn)
J1 = (10^-9)*[9060235 0 0;
              0 9060235 0;
              0 0 3654338];  %kg*m^2 main body
J2 = J1/2;                  %end body
m1 = 12;                    %kg main body
m2 = 4;                     %kg end body
M = m1+m2;                  %combined mass
mu = m1*m2/M;               %mu = m1m2/M for simplification

% Initial conditions (arbitrarily set for now)
l0 = 5;                     %initial tether length (m)
ldot0 = 0.15;                %initial tether release speed (m/s)
wx0 = 0.03;                    %initial angular velocities (rad/s)
wy0 = 0.05;
wz0 = 0.0;
l_stop = 30;                 %length of tether at end of simulation (m)


%set up state space
X = [];
X(1:3,1) = [wx0;wy0;wz0];
X(4,1) = l0;
X(5:7,1) = [0;0;0];         %omega dot initially 0
X(8,1) = ldot0;             %initial tether release speed
X(9:11,1) = [0;0;0];        %3-2-1 Euler angles initially Nadir-pointing
X(12:14,1) = [0,0,0];       %initially steady flight

%% Dynamic Model
f = @(J,w,dJdl,ldot,Mc) (J\(Mc-til(w)*J*w - ldot*dJdl*w));
%J is the inertia matrix of the tethered system about its CG
%w is the angular rate of the tethered system in body x,y,z frame
%dJdl is the rate of change of the inertia of the system with respect to tether length
%ldot is the rate of tether release (negative is contraction)
%Mc is the net moment about applied to the system about the CG (body x,y,z)
g = @(wBody,phi,theta) (inv_C(phi,theta)*wBody);
%additionally, we model the tether deployment speed as constant. In
%reality, we are able to control this in conjunction with control torques
%to produce a better response.

%% RK4 or FE
i = 1;
t = [];
t(1) = 0;

if l_stop <= l0 || ldot0 <= 0
    error('Cannnot simulate tether retraction');
end

while X(4,i) < l_stop 
    F = [0,0,0];
    r = [0,0,0];
    controller_inputs = [];
    controller_params = [];
    M_ctl = calcControlTorque(controller_params,controller_inputs);
    M_other = [0,0,0];
    pointMoments = [M_ctl; M_other];
    Mc = calcMc(F,r,pointMoments);
    if method == 1
        %FE
        ldotk1 = X(8,i);                        %current tether release speed. To be calculated by controller
        wk1 = X(1:3,i);
        lk1 = X(4,i);                           %current tether length
        r1k1 = m2/M*lk1; r2k1 = m1/M*lk1;       %current distances of end bodies from CG   
        %current inertia matrix of spacecraft
        Jk1 = J1 + m1*[r1k1^2 0 0; 0 r1k1^2 0; 0 0 0] + J2 + m2*[r2k1^2 0 0; 0 r2k1^2 0; 0 0 0]; 
        dJdlk1 = [2*lk1*mu 0 0; 0 2*lk1*mu 0; 0 0 0];    %current change in inertia wrt tether length
        k1 = dt*f(Jk1,wk1,dJdlk1,ldotk1,Mc);
        k1_l = dt*ldotk1;
        X(1:3,i+1) = wk1 + k1;
        X(4,i+1) = lk1 + k1_l;
        X(5:7,i+1) = f(Jk1,wk1,dJdlk1,ldotk1,Mc);
        X(8,i+1) = ldotk1;
        
        phik1 = X(9,i);
        thetak1 = X(10,i);
        psik1 = X(11,i);
        wBody = [wk1(1); wk1(2); wk1(3)];
        EAdotk1 = g(wBody,phik1,thetak1);
        
        
        X(9:11,i+1) = X(9:11,i) + EAdotk1 * dt;
        X(12:14,i+1) = EAdotk1;
    else
    
        %RK4

        %k1:
        ldotk1 = X(8,i);                        %current tether release speed. To be calculated by controller
        wk1 = X(1:3,i);
        lk1 = X(4,i);                           %current tether length
        r1k1 = m2/M*lk1; r2k1 = m1/M*lk1;       %current distances of end bodies from CG   
        %current inertia matrix of spacecraft
        Jk1 = J1 + m1*[r1k1^2 0 0; 0 r1k1^2 0; 0 0 0] + J2 + m2*[r2k1^2 0 0; 0 r2k1^2 0; 0 0 0]; 
        dJdlk1 = [2*lk1*mu 0 0; 0 2*lk1*mu 0; 0 0 0];    %current change in inertia wrt tether length
        k1 = dt*f(Jk1,wk1,dJdlk1,ldotk1,Mc);
        k1_l = dt*ldotk1;

        %k2:
        ldotk2 = X(8,i);                       
        wk2 = wk1 + k1/2;
        lk2 = lk1 + k1_l/2;
        r1k2 =  m2/M*lk2; r2k2 = m1/M*lk2;
        Jk2 = J1 + m1*[r1k2^2 0 0; 0 r1k2^2 0; 0 0 0] + J2 + m2*[r2k2^2 0 0; 0 r2k2^2 0; 0 0 0];
        dJdlk2 = [2*lk2*mu 0 0; 0 2*lk2*mu 0; 0 0 0];
        k2 = dt*f(Jk2,wk2,dJdlk2,ldotk2,Mc);
        k2_l = dt*ldotk2;

        %k3:
        ldotk3 = X(8,i);                       
        wk3 = wk1 + k2/2;
        lk3 = lk1 + k2_l/2;
        r1k3 =  m2/M*lk3; r2k3 = m1/M*lk3;
        Jk3 = J1 + m1*[r1k3^2 0 0; 0 r1k3^2 0; 0 0 0] + J2 + m2*[r2k3^2 0 0; 0 r2k3^2 0; 0 0 0];
        dJdlk3 = [2*lk3*mu 0 0; 0 2*lk3*mu 0; 0 0 0];
        k3 = dt*f(Jk3,wk3,dJdlk3,ldotk3,Mc);
        k3_l = dt*ldotk3;

        %k4:
        ldotk4 = X(8,i);
        wk4 = wk1 + k3;
        lk4 = lk1 + k3_l;
        r1k4 =  m2/M*lk4; r2k4 = m1/M*lk4;
        Jk4 = J1 + m1*[r1k4^2 0 0; 0 r1k4^2 0; 0 0 0] + J2 + m2*[r2k4^2 0 0; 0 r2k4^2 0; 0 0 0];
        dJdlk4 = [2*lk4*mu 0 0; 0 2*lk4*mu 0; 0 0 0];
        k4 = dt*f(Jk4,wk4,dJdlk4,ldotk4,Mc);
        k4_l = dt*ldotk4;

        %final prediction
        ldotnext = ldotk4;
        wnext = wk1 + 1/6*(k1 + 2*k2 + 2*k3 + k4);
        lnext = lk1 + 1/6*(k1_l + 2*k2_l + 2*k3_l + k4_l);
        r1next =  m2/M*lnext; r2next = m1/M*lnext;
        Jnext = J1 + m1*[r1next^2 0 0; 0 r1next^2 0; 0 0 0] + J2 + m2*[r2next^2 0 0; 0 r2next^2 0; 0 0 0];
        dJdlnext = [2*lnext*mu 0 0; 0 2*lnext*mu 0; 0 0 0];
        wdotnext = dt*f(Jnext,wnext,dJdlnext,ldotnext,Mc);


        X(1:3,i+1) = wnext;
        X(4,i+1) = lnext;
        X(5:7,i+1) = wdotnext;
        X(8,i+1) = ldotnext;
    end
    
    t(i+1) = t(i) + dt;
    i = i+1;
    
    
end
    
    
%% Plotting
figure(1); hold on;
plot(X(4,:),X(1,:));
plot(X(4,:),X(2,:));
plot(X(4,:),X(3,:));
legend('wx', 'wy', 'wz');
title('Angular Velocities in Body Fame');
xlabel('Tether Length (m)');
ylabel('Angular Velocity (rad/s)');
grid on;

figure(2); hold on;
plot(X(4,:),X(5,:));
plot(X(4,:),X(6,:));
plot(X(4,:),X(7,:));
legend('\alphax', '\alphay', '\alphaz');
title('Angular Accelerations in Body Fame');
xlabel('Tether Length (m)');
ylabel('Angular Acceleration (rad/s^2)');

figure(3); hold on;
plot(X(4,:),X(9,:));
plot(X(4,:),X(10,:));
plot(X(4,:),X(11,:));
legend('\phi body', '\theta body', '\psi body');
title('Body Frame Angles');
xlabel('Tether Length (m)');
ylabel('Anglular Displacement (rad)');

%{
figure(4); hold on;
%This is wrong
xangle_ground = X(9,:);
yangle_ground = X(10,:);
zangle_ground = X(11,:);
plot(X(4,:),xangle_ground);
plot(X(4,:),yangle_ground);
plot(X(4,:),zangle_ground);
legend('\angle x', '\angle y', '\angle z');
title('Angular Position Relative to Inertial Orbit Frame');
xlabel('Tether Length (m)');
ylabel('Anglular Position (rad)');
%}

clear; close all; clc;
global dt;
global G;
global mu;
global m_sat;
global m_pico;
global m_e;
global r1_orbit;                %initial orbital height of main body
global r2_orbit;                %initial orbital height of picosat
global w1;                      %initial orbital angular velocity of main body
global w2;                      %initial orbital angular velocity of picosat
global T_orbit;
global E;
global tether_A;
dt = 0.03;
%spacecraft and environment specs
m_sat = 1;                      %kg main body mass
m_pico = 0.1;                   %kg picosat mass
G = 6.67408*10^-11;             %m^3 kg^-1 s^-2 gravitational constant
m_e = 5.972*10^24;              %kg mass of earth
mu = G*m_e;                     %m^3 s^-2 reduced mass
E = 124 * 10^9;                 %Tether modulus (Pa)
tether_A = pi*(1*10^-4/2)^2;    %Tether cross sectional area (m^2)
l_0 = 10;                        %initial tether length (assumed deployed in the nadir direction)

%orbit specs
r1_orbit = (6378 + 350)*1000;           %initial orbit radius of main body
w1 = [0; -sqrt(mu/(r1_orbit^3)); 0];    %orbit angular velocity (rad/s)
T_orbit = 2*pi/-w1(2);                  %orbit period (s)
r2_orbit = r1_orbit - l_0;              %initial orbital radius of end body
w2 = [0; -sqrt(mu/(r2_orbit^3)); 0];    %orbital angular velocity

%Initial conditions in the orbit frame
r1 = [0;0;r1_orbit];               %initially same in inertial frame
r2 = [0;0;r2_orbit];                  %subtract initial tether length
T1 = 0.2 * r1./norm(r1);             %initial spring force (N)
T2 = -T1;
r1dot = [w1(2)*r1_orbit;0;0];                  %different in inertial frame
r2dot = [w2(2)*r2_orbit;0;0];
r1ddot = [0;0;0];
r2ddot = [0;0;0];
ldot = norm(T1)*dt/m_pico;          %impulse imparted on picosat converts to tether kickoff velocity

%convert initial conditions to inertial frame where needed
r1dotN = r1dot + cross(w,r1);
r2dotN = r2dot + cross(w,r2);
r1ddotN = r1ddot + 2.*cross(w,r1dot) + cross(w,cross(w,r1));    %wdot = 0
r2ddotN = r2ddot + 2.*cross(w,r2dot) + cross(w,cross(w,r2));


t = 0;
k = 1;
l = l_0;
while t < 20
    %compute next state
    [r1k1, r2k1, r1dotk1, r2dotk1, r1ddotk1, r2ddotk1] = satmotion(r1(:,k),r2(:,k),r1dotN(:,k),r2dotN(:,k),T1(:,k),T2(:,k));
    
    %compute tether release length
    lk1 = ldot*dt + l(k); % assume constant tether release speed                
    
    %compute tension
    
    rhatk1 = r2k1-r1k1 ./ norm(r2k1-r1k1);      %unit direction vector for tension
    r12N = r2k1 - r1k1;
    
    T1k1 = E*tether_A*(norm(r12N) - l(k))./lk1 .* rhatk1; % get diff in current pos and tether length returns a strain that
    T2k1 = -T1k1;                                           % we can use to compute the tension


    %update state
    r1(:,k+1) = r1k1;
    r2(:,k+1) = r2k1;
    r1dotN(:,k+1) = r1dotk1;
    r2dotN(:,k+1) = r2dotk1;
    r1ddotN(:,k+1) = r1ddotk1;
    r2ddotN(:,k+1) = r2ddotk1;
    T1(:,k+1) = T1k1;
    T2(:,k+1) = T2k1;
    l(k+1) = lk1;
    
    
    t = t+dt;
    k = k+1;
end







%% Satellite Dynamics Model
%Simple model to describe motion of each end body satellite in orbit. Each
% end body is modeled as a point mass, subjected to the gravitational pull
% of the earth and another net force (denoted T1, T2). To solve for the
% motion of each body, we integrate the two body problem with the
% additional T1 T2 force.
function [r1_out,r2_out,r1dot_out,r2dot_out,r1ddot,r2ddot] = satmotion(r1,r2,r1dot,r2dot,T1,T2)
%state equation describing motion of each satellite
%basically F = ma marching
%inputs:
%   r1, r2 are the position vectors of each body in the inertial frame
%   r1dot, r2dot are the first time derivatives of the position vectors of
%                each body in the  inertial frame
%   T1, T2 are the force vectors due to all forces except gravity for each
%          body (also in inertial frame)
%outputs:
%   r1_out, r2_out are the r1 and r2 of the next state (inertial frame)
%   r1dot_out, r2dot_out are r1dot and r2dot of the next state (inertial f)
%   r1ddot, r2ddot are the second time derivatives of position vectors of
%                  each body (inertial frame)
global mu;
global dt;
global m_sat;
global m_pico;
%unit direction vectors in inertial frame
r1hat = r1./norm(r1);
r2hat = r2./norm(r2);

%F = ma for each body
r1ddot = -mu./(norm(r1)^2).*r1hat;% + T1./m_sat;    %Testing without T
r2ddot = -mu./(norm(r2)^2).*r2hat;% + T2./m_pico;

%integrating rddot to rdot
r1dot_out = r1dot + r1ddot.*dt;
r2dot_out = r2dot + r2ddot.*dt;

%integrating rdot to r
r1_out = r1 + r1dot_out.*dt;
r2_out = r2 + r2dot_out.*dt;
end
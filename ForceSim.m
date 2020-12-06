clear; close all;
%Coordinate definition: xyz orthoganal, inertial. Orbit plane is x-y plane
%Assumed that deployment begins when spacecraft is on the positive y-axis,
%   and the spacecraft moves counterclockwise in the x-y plane
global dt;
global G;
global Re;
global mu;
global m_sat;
global m_pico;
global m_e;
global E;
global tether_A;
dt = 0.2;
%known parameters
m_sat = 1;                      %kg main body mass
m_pico = 0.270;                   %kg picosat mass
Re = 6378e3;                    %m earth radius
G = 6.67408e-11;                %m^3 kg^-1 s^-2 gravitational constant
m_e = 5.972e24;                 %kg mass of earth
mu = G*m_e;                     %m^3 s^-2 reduced mass
E = 117e9;                      %Tether modulus (Pa)
tether_A = pi*(1e-3/2)^2;    %Tether cross sectional area (m^2)
l_0 = 0.1;                        %initial tether length (assumed deployed in the nadir direction)
ldot = 1;                       %m/s tether release rate
%Inital conditions
r1_init = [0; Re+350e3; 0];     %Main body orbit at 350 km altitude
r2_init = [0; Re+350e3-l_0; 0];     %picosat orbit at at 350 km altitude (attached to main body)
w1_init = [0; 0; sqrt(mu/norm(r1_init)^3)];
w2_init = [0; 0; sqrt(mu/norm(r2_init)^3)];
r1dot_init = cross(w1_init, r1_init);
r2dot_init = [0; -ldot; 0] + cross(w2_init, r2_init);
r1ddot_init = -mu.*r1_init./norm(r1_init)^3;
r2ddot_init = -mu.*r2_init./norm(r2_init)^3;

T_init = [0; 0; 0];

t = 0;
k = 1;
l = l_0;
r1 = r1_init; r2 = r2_init;
r1dot = r1dot_init; r2dot = r2dot_init;
r1ddot = r1ddot_init; r2ddot = r2ddot_init;
T = T_init;


while t < 2*pi/norm(w1_init)
    %predict next state
    %[r1_pred,r2_pred] = satmotion(r1(:,k),r2(:,k),r1dot(:,k),r2dot(:,k),[0;0;0]);

    %compute tether release length
    lk1 = ldot*dt + l(k);  
                   
    
    r12 = r1(:,k)-r2(:,k);
    r12hatk1 = r12./norm(r12);
    delta = norm(r12) - l(k);
   
    %tether carries no compressive strain (slack tether case)
    if delta < 0
        delta = 0;
    end
   
    
    T(:,k) = E*tether_A*delta./l(k) .* r12hatk1;

    %compute next state
    [r1k1, r2k1, r1dotk1, r2dotk1] = satmotion(r1(:,k),r2(:,k),r1dot(:,k),r2dot(:,k),T(:,k));
    
    %update state
    r1(:,k+1) = r1k1;
    r2(:,k+1) = r2k1;
    r1dot(:,k+1) = r1dotk1;
    r2dot(:,k+1) = r2dotk1;
    l(k+1) = lk1;
    
    
    t = t+dt;
    k = k+1;
end
plot(r1(1,:),r1(2,:)); pbaspect([1 1 1]); hold on;
plot(r2(1,:),r2(2,:)); pbaspect([1 1 1]);



%% Satellite Dynamics Model
%Simple model to describe motion of each end body satellite in orbit. Each
% end body is modeled as a point mass, subjected to the gravitational pull
% of the earth and another net force (denoted T1, T2). To solve for the
% motion of each body, we integrate the two body problem with the
% additional T1 T2 force.
function [r1_out,r2_out,r1dot_out,r2dot_out] = satmotion(r1,r2,r1dot,r2dot,T)
%state equation describing motion of each satellite
%
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
r1ddot = -mu./(norm(r1)^2).*r1hat - T./m_sat;
r2ddot = -mu./(norm(r2)^2).*r2hat + T./m_pico;

%integrating rddot to rdot
r1dot_out = r1dot + r1ddot.*dt;
r2dot_out = r2dot + r2ddot.*dt;

%integrating rdot to r
r1_out = r1 + r1dot_out.*dt;
r2_out = r2 + r2dot_out.*dt;
end

%% Two body prediction with no tether
% currently not in use

function[r1_pred, r2_pred] = predict(r1,r2,ldot)
global mu;
global dt;
w1 = [0; 0; sqrt(mu/norm(r1)^3)];
w2 = [0; 0; sqrt(mu/norm(r2)^3)];

r1dot = cross(w1, r1);
r2dot = -ldot.*r2./norm(r2) + cross(w2, r2);

r1_pred = r1 + dt.*r1dot;
r2_pred = r2 + dt.*r2dot;
end
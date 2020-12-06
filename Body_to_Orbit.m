function [orbit_angles] = Body_to_Orbit(body_angles)

a = body_angles(1); %roll
b = body_angles(2); %pitch
c = body_angles(3); %yaw
DCM = [cos(b)*cos(c) cos(b)*sin(c) -sin(b)


end
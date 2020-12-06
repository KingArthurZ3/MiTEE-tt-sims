function[Mc] = calcMc(F, r, M)
%calcMc calculates the net moment about the center of gravity of an object
%Inputs:
%   F - nx3 matrix who's columns are the x,y,z components of each force
%       (each row is a force)
%   r - nx3 matrix who's columns are the x,y,z components of each force's
%       point of application (with respect to CG, in the body frame)
%   M - mx3 matrix who's columns are the x,y,z components of each point
%       moment applied about the CG
%Outputs:
%   Mc - 3x1 vector of the net moment about the CG in body x,y,z frame
%All forces,  moments, and positions should be in body x,y,z coordinates

F = F'; r = r'; M = M';
all_moments = [cross(r,F),M];
Mc = sum(all_moments,2);

end
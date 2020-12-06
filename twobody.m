function Xdot = twobody(t,X)
%   X: [rx; ry; rz; rx'; ry'; rz']      
%   Xdot: [rx'; ry'; rz'; rx''; ry''; rz'']
%in inertial coordinate frame

global mu;
r_mag = norm(X(1:3));
Xdot(1:3) = X(4:6);
Xdot = Xdot';
Xdot(4:6) = -mu/(r_mag^3) .* X(1:3);

end
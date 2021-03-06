function invc = inv_C(phi,theta)
%inverse cosine matrix
invc = 1/cos(theta)*[cos(theta) sin(phi)*sin(theta) cos(phi)*sin(theta);
                        0 cos(phi)*cos(theta) -sin(phi)*cos(theta);
                        0 sin(phi) cos(phi)];
end


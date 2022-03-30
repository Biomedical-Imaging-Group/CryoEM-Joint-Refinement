function [phi, theta, psi] = q_to_euler(q)
%% Converts the quaternion representation of rotation to euler angles in zyz
% convention
% used the formulations in this web link:...
% http://bediyap.com/programming/convert-quaternion-to-euler-rotations/
% Mona Zehni, 2020

phi = atan2( 2 *(q(3)*q(4)-q(1)*q(2)), 2*(q(2)*q(4)+q(1)*q(3)) );
theta = acos ( q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2 );
psi = atan2( 2*(q(3)*q(4) + q(1)*q(2)), -2*(q(2)*q(4)-q(1)*q(3)) );

end
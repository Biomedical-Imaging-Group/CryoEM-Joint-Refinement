%=========================================================================%
% This function randomly generates equidistributed proj angles. 
% Psi is automatically set at zero. 
% Author: laurene.donati@epfl.ch
% BASED ON DESERNO PAPER. 
% Self-note: In our notations, his theta is our tilt, his phi is our rot
%-------------------------------------------------------------------------%

function [rot, tilt, psi] = generateEquidistribRandomProjAngles(nberProj)

 % Set rot 
 z = (rand(1,nberProj)*2)-1;
 rot  = rand(1,nberProj)*2*pi;
 
 % Set tilt 
 x    = sqrt(1-(z.^2)).*cos(rot);
 y    = sqrt(1-(z.^2)).*sin(rot);
 tilt = atan(sqrt(y.^2+x.^2)./z);
 tilt(tilt<0) = tilt(tilt<0)+pi; 
 
 % Set psi 
 psi = rand(1, nberProj) * 2 * pi;
 
end 

%%    Copyright (C)
%     2018 - Biomedical Imaging Group, EPFL
%     Masih Nilchian - masih.nilchian@epfl.ch
%     Laurène Donati - laurene.donati@epfl.ch
%     Slightly modified by Mona Zehni

function [r1,r2] = setProjectionPlaneOrientationsSamples(angles, rotMatrix)
num_samples = size(angles,1);
rot = angles(:,1);
tilt = angles(:,2);
psi = angles(:,3);
sizeAngles = num_samples;
r1 = zeros(3,sizeAngles);
r2 = zeros(3,sizeAngles);

for n = 1:num_samples
    R = Euler(rot(n), tilt(n), psi(n));
    R = R * rotMatrix;
    r1(:,n) = [R(1,1),R(1,2),R(1,3)]';
    r2(:,n) = [R(2,1),R(2,2),R(2,3)]';
end

end

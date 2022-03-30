%%    Copyright (C) 
%     2018 - Biomedical Imaging Group, EPFL
%     Masih Nilchian - masih.nilchian@epfl.ch	
%     Laurène Donati - laurene.donati@epfl.ch


function [r1,r2] = setProjectionPlaneOrientations(rot,tilt,psi)

sizeAngles = length(rot)*length(tilt)*length(psi); 
r1 = zeros(3,sizeAngles); 
r2 = zeros(3,sizeAngles); 
count  = 1;
for i = 1 :length(rot)
    for j = 1 : length(tilt)
        for k = 1 : length(psi)
            % Account for different addressing conventions
            R = Euler(rot(i), tilt(j), psi(k));
            r1(:,count) = [R(1,1),R(1,2),R(1,3)]';
            r2(:,count) = [R(2,1),R(2,2),R(2,3)]';
            count = count+1;
        end
    end
end

end

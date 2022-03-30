function [rotNew, tiltNew, psiNew, rot, tilt, psi] = align_relion_angles(volRelion, data_star, numProj, gt_angles, volGT, rotMatrix)
%% Align the angles output from Relion with the convention used in LinOpPBTShift
% check the consistencies between the projection view conventions

%input volRelion: 3D ab inito volume from Relion
%input data_star: path to the relion star file containig the intial
%projection and in-plane translations estimated from Relion
%input numProj: # of projections
%input gt_angles: ground truth (GT) angles used to generate the projections usin
%g LinOpPBTShift
%input volGT: GT volume
%input rotMatrix: 3x3 rotation matrix used to align volRelion and volGT

%output: aligned projection angles

% read the anlges
data = readSTAR(data_star);
for t = 1:numProj
    t
    rotFull(t) = data.data{t}.rlnAngleRot; 
    tiltFull(t) = data.data{t}.rlnAngleTilt; 
    psiFull(t) = data.data{t}.rlnAnglePsi;
end

ind = 200;
rot =  ((90 - rotFull)/360)*2*pi;
tilt = (tiltFull/360)*2*pi;
psi = (-(90 + psiFull)/360)*2*pi;

H = LinOpPBTShift(size(volGT), [gt_angles(ind, 1), gt_angles(ind, 2), gt_angles(ind, 3)], [0, 0], 1, 1, 1, 1);
projGT = H.apply(volGT);

H = LinOpPBTShift(size(volRelion), [rot(ind), tilt(ind), psi(ind)], [0, 0], 1, 1, eye(3), 1, 1);
projRelion = H.apply(volRelion);

figure;
subplot(1,2,1); imagesc(projRelion); title('Relion');
subplot(1,2,2); imagesc(projGT); title('LinOpPBTShift');

for i = 1:numProj
    q = rot_to_q((Euler(rot(i), tilt(i), psi(i))*(rotMatrix)).');
%     q(2) = -q(2); q(3) = -q(3);
    [rotNew(i), tiltNew(i), psiNew(i)] = q_to_euler(q);
end

H = LinOpPBTShift(size(volGT), [rotNew(ind), tiltNew(ind), psiNew(ind)], [0, 0], 1, 1, eye(3), 1, 1);
projPBT1 = H.apply(volGT);

H = LinOpPBTShift(size(volGT), [gt_angles(ind, 1), gt_angles(ind, 2), gt_angles(ind, 3)], [0, 0], 1, 1, eye(3), 1, 1);
projPBTOrig = H.apply(volGT);

figure;
subplot(1, 3, 1); imagesc(projRelion); title('Proj Relion');
subplot(1, 3, 2); imagesc(projPBTOrig); title('PBTOrig');
subplot(1, 3, 3); imagesc(projPBT1); title('PBT1');

end
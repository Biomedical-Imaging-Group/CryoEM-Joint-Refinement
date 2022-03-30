function [estR, volAligned, reflect] = align_dl_densities(volRef, vol)
%% Align vol with volRef and outputs the rotation matrix and the aligned
%volume.
%input volRef: Reference volume
%input vol: volume to be aligned to ref

%output estR: rotation matrix that aligns vol with volRef
%output volAligned: aligned volume
%output reflect: whether there is a reflection needed for alignment

% Mona Zehni, March 2022

% to allow faster alignment, one can downsample the volumes
dl = 6;
volRefDL = volRef(1:dl:size(volRef, 1), 1:dl:size(volRef, 2), 1:dl:size(volRef, 3));
volDL = vol(1:dl:size(vol, 1), 1:dl:size(vol, 2), 1:dl:size(vol, 3));
[estR, ~, volAligned, reflect] = cryo_align_densities(volRefDL, volDL);

end


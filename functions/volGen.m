function [vol, sz] = volGen(filename, DL)
%% Downsamples the volume based on the DL rate
%input filename: name of the file containing the volume
%input DL: the downsample factor (integer)
%output: the downsampled volume
% Mona Zehni, July 2018

vol = double(ReadMRC(filename));

% if the DL  > 1, consider low-pass filtering of vol before down-sampling

% manually downsampling the volume
rows = 1:DL:size(vol,1);
cols = 1:DL:size(vol,2);
zs = 1:DL:size(vol,3);

vol = vol(rows,:,:);
vol = vol(:,cols,:);
vol = vol(:,:,zs);
sz = size(vol);
end


function vol_dl = downSample(vol, scaleRatio)
%% Downsamples the volume
%input vol: the volume
%input scaleRatio: downsample factor
%output vol_dl: Downsampled volume

% Mona Zehni, July 2018

% if the scaleratio  > 1, consider low-pass filtering before calling
% downSample
rows = scaleRatio:scaleRatio:size(vol,1);
cols = scaleRatio:scaleRatio:size(vol,2);
zs = scaleRatio:scaleRatio:size(vol,3);

vol_dl = vol;
vol_dl = vol_dl(rows,:,:);
vol_dl = vol_dl(:,cols,:);
vol_dl = vol_dl(:,:,zs);

end
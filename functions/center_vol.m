function [vol_pad, cx, cy, cz] = center_vol(vol_sz, vol)
%% Centralizes the volume
%param vol_sz: final volume size
%param vol: density map from Chimera (grid spacing 0.5)

%return vol_pad: centralized zero padded volume
%return cx, cy, cz: final center of masses, should be close to [0, 0, 0]
% Mona Zehni, 2021

n = floor(vol_sz/2);
c = ceil(vol_sz/2);
[N1, N2, N3] = size(vol);
index_x = find_index(c, N1);
index_y = find_index(c, N2);
index_z = find_index(c, N3);

vol_pad = zeros(vol_sz, vol_sz, vol_sz);
vol_pad(index_x, index_y, index_z) = vol;
[cx, cy, cz] = compute_cent_mass(n, vol_pad);

vol_pad = zeros(vol_sz, vol_sz, vol_sz);
vol_pad(index_x-cx, index_y-cy, index_z-cz) = vol;
% the center of mass should be close to [0, 0, 0]
[cx, cy, cz] = compute_cent_mass(n, vol_pad);
end

function [cx, cy, cz] = compute_cent_mass(n, vol)
if mod(size(vol, 1), 2) == 1
    [X, Y, Z] = meshgrid(-n:n, -n:n, -n:n);
else
    [X, Y, Z] = meshgrid(-n:n - 1, -n:n - 1, -n:n - 1);
end
cx = vol.*X;
cx = ceil(sum(cx(:))/sum(vol(:)));
cy = vol.*Y;
cy = ceil(sum(cy(:))/sum(vol(:)));
cz = vol.*Z;
cz = ceil(sum(cz(:))/sum(vol(:)));
end

function index = find_index(c, N)
if mod(N, 2)==1
    index = c-floor(N/2):c+floor(N/2);
else
    index = c-floor(N/2):c+floor(N/2)-1;
end
end
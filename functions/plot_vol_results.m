function plot_vol_results(vol1, vol2, title1, title2)
%% Plots vol1 and vol2 with the corresponding titles
% Mona Zehni, March 2022

vol_sz = size(vol1, 1);

figure;
subplot(2, 2, 1)
view3d(vol1)
title(title1)
subplot(2, 2, 2)
view3d(vol2)
title(title2)
subplot(2, 2, 3)
imagesc(vol1(:, :, floor(vol_sz/2)));
colorbar
ylabel('Central slice')
colormap gray
subplot(2, 2, 4)
imagesc(vol2(:, :, floor(vol_sz/2)))
colorbar
colormap gray

end
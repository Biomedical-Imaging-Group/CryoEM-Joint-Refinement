function plot_hist(gt_angles, gt_shifts, est_angles, est_shifts, angles_init, shifts_init)
%% Plots the hist of the error in the estimated angles/shifts at different refinement iterations

%input gt_angles: ground truth (GT) angles
%input gt_shifts: ground truth (GT) shifts
%input est_angles: estimated angles for different iterations
%input est_shifts: estimated shifts for different iterations
%input angles_init: initial angles
%input shifts_init: initial shifts

num_iter = size(est_angles, 3);

figure;
% plot angles
for i = 0:6:num_iter
    if i == 0
        tmp = angles_init - gt_angles;
    else
        tmp = est_angles(:, :, i) - gt_angles;
    end
    for j = 1:1
        
        tmp_ang = tmp(:, j);
        tmp_ang = min(tmp_ang, 2 * pi - tmp_ang);
        tmp1 = 2 * pi + tmp_ang;
        tmp_ang = tmp_ang .* (abs(tmp_ang) < abs(tmp1)) + tmp1 .* (abs(tmp1) < abs(tmp_ang));
        
        [N, X] = hist(tmp_ang, 30);
        sum(N)
        plot(X, N/sum(N))
        hold on
        
    end
end

% plot shifts
figure
for i = 0:6:num_iter
    if i == 0
        tmp = shifts_init - gt_shifts;
    else
        tmp = est_shifts(:, :, i) - gt_shifts;
    end
    for j = 1:1
        tmp = tmp(:, j);
        
        [N, X] = hist(tmp, 30);
        sum(N)
        plot(X, N/sum(N)) 
        hold on
        
    end
end



end
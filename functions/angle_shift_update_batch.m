function [updated_angle, updated_shift] = ...::
    angle_shift_update_batch(vol, proj, angles_init, shifts_init, angle_step, ...
    shift_step, learning_rate_angle, learning_rate_shift, alpha, beta, DL_proj)
%Updating the unknown parameters of angles and shifts based on gradient
%steps plus using line search method

% Note that this file compared to the back-up file take into account the
% in-plane rotations while the bu file works for the case where we have no
% in-plane rotations.
%Input:
%       vol: the current reconstructed volume
%       proj: the set of projection images
%       angles_init: the initial angles
%       shifts_init: the initial shifts
%       angle_step: the angle step used in order to compute the gradients
%       wrt rotations and tilts numerically
%       shift_step: the shift step used in order to compute the gradients
%       wrt shifts x and y numerically
%       learning_rate_angle: the starting learning rate for the angles
%       learning_rate_shift: the starting learning rate for the shifts
%       alpha, beta: the parameters of the line search
%Output:
%       updated_angle: the updated angles after performing one step of
%       gradient descent with line search
%       updated_shift: the updated shifts after performing one step of
%       gradient descent with line search
% Mona Zehni, July 2018

num_samples = size(proj,3);

proj_only = 1;
scaleVal = 1;

kbwf_config.a = 4;
kbwf_config.alpha = 19;
kbwf_config.m = 2;
kbwf_config.scale = 1;
kbwf_config.Ti = 1;
kbwf_config.Tp = 1;

% limiting the while loop for the line search, for speed up purposes.
iter_max_ls = 3; % it was 5 previously

% initializing the variables
rot_updated = angles_init(:,1);
tilt_updated = angles_init(:,2);
in_plane_updated = angles_init(:,3);

shiftX_updated = shifts_init(:,1);
shiftY_updated = shifts_init(:,2);

% gradient wrt the rotations
H_init = LinOpPBTShift(size(vol),angles_init,shifts_init,proj_only,kbwf_config);
proj_init = H_init.apply(vol);
cost1 = reshape(sum(sum((proj-proj_init).^2,1),2),[num_samples,1]);

H_rot = LinOpPBTShift(size(vol),[angles_init(:,1)+angle_step,angles_init(:,2:3)],shifts_init,proj_only,kbwf_config);
proj_rot = H_rot.apply(vol);
diff_rot = 2*sum(sum((proj_rot-proj_init).*(proj_init-proj),1),2);
deriv_rot = reshape(diff_rot / angle_step,[num_samples,1]);

% gradient wrt tilts
H_tilt = LinOpPBTShift(size(vol),[angles_init(:,1),angles_init(:,2)+angle_step,angles_init(:,3)],shifts_init,proj_only,kbwf_config);
proj_tilt = H_tilt.apply(vol);
diff_tilt = reshape(2*sum(sum((proj_tilt-proj_init).*(proj_init-proj),1),2),[num_samples,1]);
deriv_tilt = reshape(diff_tilt / angle_step, [num_samples,1]);

% gradient wrt in-plane rotations
H_in_plane = LinOpPBTShift(size(vol),[angles_init(:,1),angles_init(:,2),angles_init(:,3)+angle_step],shifts_init,proj_only,kbwf_config);
proj_in_plane = H_in_plane.apply(vol);
diff_in_plane = reshape(2*sum(sum((proj_in_plane-proj_init).*(proj_init-proj),1),2),[num_samples,1]);
deriv_in_plane = reshape(diff_in_plane / angle_step, [num_samples,1]);

% updating the angles
lr = learning_rate_angle;
ind_acc = []; % the set containnig the indices of the updated projection angles
index = [1:num_samples];
rot_tmp = angles_init(:,1);
tilt_tmp = angles_init(:,2);
in_plane_tmp = angles_init(:,3);
all_updated_rot = 0;
iter = 1;
while all_updated_rot==0
    % applying the gradients over the angles
    H_angle = LinOpPBTShift(size(vol),...
        [rot_tmp(index)-lr*deriv_rot(index),tilt_tmp(index)-lr*deriv_tilt(index), in_plane_tmp(index)-lr*deriv_in_plane(index)],shifts_init(index,:),proj_only,kbwf_config);
    
    proj_angle = H_angle.apply(vol);
    cost2 = reshape(sum(sum((proj(:,:,index)-proj_angle).^2,1),2),[length(index),1]);
    
    ind = find(((cost2-cost1(index))<0)); 
    rot_updated(index(ind),1) = mod(rot_tmp(index(ind)) - lr*deriv_rot(index(ind)),2*pi);
    tilt_updated(index(ind),1) = mod(tilt_tmp(index(ind)) - lr*deriv_tilt(index(ind)),2*pi);
    in_plane_updated(index(ind),1) = mod(in_plane_tmp(index(ind)) - lr*deriv_in_plane(index(ind)),2*pi);
    index(ind) = [];
    ind_acc = [ind_acc;ind];
    
    all_updated_rot = (length(ind_acc)==num_samples) | (iter>iter_max_ls);
    iter = iter + 1;
    lr = lr * beta;
    
end

% gradient with respect to the shifts x
angles = [rot_updated,tilt_updated,in_plane_updated];
H_init = LinOpPBTShift(size(vol),angles,shifts_init,proj_only,kbwf_config);
proj_init = H_init.apply(vol);
cost1 = reshape(sum(sum((proj-proj_init).^2,1),2),[num_samples,1]);

H = LinOpPBTShift(size(vol),angles,[shifts_init(:,1)+shift_step,shifts_init(:,2)],proj_only,kbwf_config);
proj_shiftX = H.apply(vol);
diff_shiftX = 2*sum(sum((proj_shiftX-proj_init).*(proj_init-proj),1),2);
deriv_shiftX = reshape(diff_shiftX / shift_step, [num_samples,1]);

H = LinOpPBTShift(size(vol),angles,[shifts_init(:,1),shifts_init(:,2)+shift_step],proj_only,kbwf_config);
proj_shiftY = H.apply(vol);
diff_shiftY = 2*sum(sum((proj_shiftY-proj_init).*(proj_init-proj),1),2);
deriv_shiftY = reshape(diff_shiftY / shift_step, [num_samples,1]);

%updating the shifts
lr = learning_rate_shift;
ind_acc = [];
index = [1:num_samples];
shiftx_tmp = shifts_init(:,1);
shifty_tmp = shifts_init(:,2);
all_updated_shiftX = 0;
iter = 1;
while all_updated_shiftX==0
    
    H_shift = LinOpPBTShift(size(vol),angles(index,:),...
        [shiftx_tmp(index)-lr*deriv_shiftX(index), shifty_tmp(index)-lr*deriv_shiftY(index)],proj_only,kbwf_config);
    proj_shift = H_shift.apply(vol);
    cost2 = reshape(sum(sum((proj(:,:,index)-proj_shift).^2,1),2),[length(index),1]);
    
    ind = find(((cost2-cost1(index))<0));
    shiftX_updated(index(ind)) = shiftx_tmp(index(ind)) - lr*deriv_shiftX(index(ind));
    shiftY_updated(index(ind)) = shifty_tmp(index(ind)) - lr*deriv_shiftY(index(ind));
    index(ind) = [];
    ind_acc = [ind_acc;ind];
    
    all_updated_shiftX = (length(ind_acc)==num_samples) | (iter>iter_max_ls) ;
    iter = iter + 1;
    lr = lr * beta;
    
end

updated_angle = [rot_updated,tilt_updated,in_plane_updated];
updated_shift = [shiftX_updated,shiftY_updated];


end

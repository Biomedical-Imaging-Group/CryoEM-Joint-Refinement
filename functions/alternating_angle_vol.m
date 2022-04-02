function [ vol, angles, shifts, final_SNR, vol_recon_cell, final_fval, final_iternum ] = ...
    alternating_angle_vol(params)
%% This function solves for the volume and the angles/shifts using alternating optimization

%input params: A struct containing the following variables:

%input max_iter: the maximum number of iterations
%input vol_init: the initial volume to start from
%input angles_init: the initial angles to start from
%input shifts_init: the initial shifts to start from
%input proj: the projection images
%input LS: the least squares operator
%input R1: the regularization operator
%input rho_n_final: the final rho_n
%input max_iter_ADMM: the maximum number of ADMM steps taken at each iteration
%input grad_num: the number of gradient decsent steps applied to update the
%angles/shifts
%input grad_step_rot: the learning rate used to update the rotation angle
%input grad_step_tilt: the learning rate used to update the tilt angle
%input alpha, beta: parameters of the line search model
%input num_scale_stages: the number of scale stages used during the
%alternating optimization, here this is always 1.
%update, here it is set as the ground truth volume.
%input lamb_init: the initial lambda used in the optimization
%input lamb_final: the final lambda used in the optimization
%input vol_mask: volume mask

%output vol: the final refined volume
%output angles, shifts: the final refined angles and shifts
%output final_SNR: the SNR of the recovered volume through iterations
%output vol_recon_cell: the cell containing the refined volumes in
%different iterations
%output final_fval: the list containing the final function values through
%iterations
%output final_iternum: the list of the iteration steps for updating the
%volume
% Mona Zehni, Aug 2018

num_pool = 2;

param_exp;

angle_step = 1e-4; % used for computing the numerical gradients
shift_step = 1e-4;
gamma = 1.7;
num_samples = size(params.projs,3);

final_SNR = [];
final_fval = [];
final_iternum = [];

% setting the stages for which we change the scales
scale_ratio = 2.^[params.num_scale_stages-1:-1:0];

num_iter = floor(1/sum(1./scale_ratio).*(1./scale_ratio)*params.max_iter);
max_iter = sum(num_iter);

% parameters of the KBWF
param = struct;
param.a = 4;
param.alpha = 19;
param.m = 2;

% decides which update, angle/volume to start from
angle_first = 1;

% grouping the projections as batches
batch_size = floor(num_samples/num_pool);
batch_ind = [1:num_samples];
batch_num = floor(num_samples/batch_size);
batch_ind = reshape(batch_ind,[batch_size,batch_num]);

% dividing the projections and the angles in batches
% we can used downsampled proj images during angle/shift updates
projs = params.projs;
proj_tmp = projs(1:params.DL_proj:size(projs,1), 1:params.DL_proj:size(projs,2), :);
proj_batch = zeros(size(proj_tmp,1),size(proj_tmp,2),batch_size,batch_num);

% only during angle/shift update we are downsampling the volume
for b = 1:batch_num
    proj_batch(:,:,:,b) = proj_tmp(:,:,batch_ind(:,b));
end

% deriving the different lambda and rho values for different steps
rho_step = (params.rho_n_final/params.rho_n_init)^(1/(params.max_iter-1));
lamb_step = (params.lamb_final/params.lamb_init)^(1/(params.max_iter-1));

rho_iter = params.rho_n_init*(rho_step.^([0:params.max_iter-1]));
lamb_iter = params.lamb_init*(lamb_step.^([0:params.max_iter-1]));

% a cell to store reconstructed volumes in different stages
vol_recon_cell = {};

% generate some of the operators
vol_init = params.vol_init;
angles_init = params.angles_init;
shifts_init = params.shifts_init;
updated_angle = zeros([size(angles_init),max_iter]);
updated_shift = zeros([size(shifts_init),max_iter]);

LS = params.LS;
R1 = params.R1;
Hn = {LinOpGrad(size(vol_init)), LinOpIdentity(size(vol_init))};
R_pos = CostNonNeg(size(vol_init));

grad_step_rot = params.grad_step_rot;
grad_step_tilt = params.grad_step_tilt;
grad_step_shift = params.grad_step_shift;

for iter = 1:params.max_iter
    
    fprintf('iter = %d\n',iter)
    rho_n = rho_iter(iter);
    lamb = lamb_iter(iter); %/1
    
    turn_angle = angle_first;
    turn_vol = ~(angle_first);
    
    % this flag shows whether a full update cycle of angle and volume is
    % done or not, when flag==2, then angles and volume are both updated.
    flag = 0;
    
    while flag < 2
        
        if turn_angle ==0 && turn_vol == 1 && flag<2
            
            % volume update
            % build the operator based on the updated angles and shifts
            fprintf('Updating volume...\n')
            H = LinOpPBTShift(size(vol_init), angles_init, shifts_init, 0, kbwf_recon);
            Fn = {lamb*R1, R_pos};
            
            [vol_recon, SNR_tmp, fval_tmp, iternum_tmp] = ...
                recon_3D_ADMM(LS, H, Fn, Hn, [rho_n/2,1e3], zeros(size(vol_init)), vol_init, params.max_iter_ADMM);
            
            % applying the mask and the positivity constraint to the
            % reconstructed volume
            vol_recon(find(vol_recon<0)) = 0;
            vol_recon = vol_recon.*params.vol_mask;
            vol_recon_cell{iter} = vol_recon;
            
            final_SNR = [final_SNR, SNR_tmp];
            final_fval = [final_fval, fval_tmp];
            % the overall number of iterations
            final_iternum = [final_iternum, (iter-1)*params.max_iter_ADMM+iternum_tmp];
            
            % update the initial volume used for the next stage
            vol_init = vol_recon;
            turn_angle = ~(turn_angle);
            turn_vol = ~(turn_vol);
            flag = flag + 1;
        end
        
        if turn_angle==1 && turn_vol==0 && flag<2
            fprintf('Updating angles/in-plane translations...\n')
            for n = 1:params.grad_num
                fprintf(['GD iter: ' num2str(n) '\n'])
                updated_angle_tmp = zeros(batch_size,3,batch_num);
                updated_shift_tmp = zeros(batch_size,2,batch_num);
                if isempty(gcp('nocreate'))
                    parpool('local',num_pool);
                end
                tic
                parfor batch_index = 1:batch_num
                    
                    [updated_angle_tmp(:,:,batch_index), updated_shift_tmp(:,:,batch_index)] = ...
                        angle_shift_update_batch...
                        (vol_init, ...
                        squeeze(proj_batch(:,:,:,batch_index)), ...
                        angles_init(batch_ind(:,batch_index),:), ...
                        shifts_init(batch_ind(:,batch_index),:), ...
                        angle_step, shift_step, grad_step_rot*10, ...
                        grad_step_shift*10, params.alpha, params.beta, ...
                        params.DL_proj);
                    
                end
                toc
                
                updated_angle_tmp = reshape(permute(updated_angle_tmp,[1,3,2]),[num_samples,3]);
                updated_shift_tmp = reshape(permute(updated_shift_tmp,[1,3,2]),[num_samples,2]);
                
                updated_angle(:,:,(iter-1)*params.grad_num+n) = updated_angle_tmp;
                updated_shift(:,:,(iter-1)*params.grad_num+n) = updated_shift_tmp;
                angles_init = updated_angle_tmp;
                shifts_init = updated_shift_tmp;
                
            end
            turn_angle = ~(turn_angle);
            turn_vol = ~(turn_vol);
            flag = flag + 1;
            
            plot_angles_shifts
        end
        %delete(gcp('nocreate'))
    end
    
    % this is done to still update the values, even when the gradient
    % becomes very small. we can remove gamma and that should lead to a
    % faster convergence.
    grad_step_rot = grad_step_rot * gamma;
    grad_step_tilt = grad_step_tilt * gamma;
    grad_step_shift = grad_step_shift * gamma;
end

vol = vol_recon;
angles = updated_angle;
shifts = updated_shift;

end

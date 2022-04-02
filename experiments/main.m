clear all
clc

rng(1);

%% loading the required parameters
param_exp;
num_samples = 200; % # of projs

%% optimization setting and definition
optim_config;

% generating the volume
% check prepare_volume function for a list of existing volumes
vol_name = 'test_3400';
vol_gt = prepare_volume(vol_name, 1);
vol_sz = size(vol_gt, 1);

fprintf('Expansion on KBWF...\n')
[vol_coeff, vol_est, kbwf_kernel] = expand_vol_KBWF(vol_gt, kbwf_proj);
plot_vol_results(vol_gt, vol_est, 'GT volume', 'Expanded volume')

% generate the angles and the shifts
sigma_shiftX = 6;
sigma_shiftY = 6;
[rot, tilt, psi] = generateEquidistribRandomProjAngles(num_samples);
shifts = ShiftGen(sigma_shiftX, sigma_shiftY, num_samples); % used to create the shifts
angles = [rot',tilt',psi'];

sigma_noise = 10.^(linspace(-1, 2, 5));
sig_noise_index = 1;
sigma_angle = 0.7; % STD of the noise on the projection angles

% forward model and projection dataset generation
H = LinOpPBTShift(size(vol_coeff), angles, shifts, 1, kbwf_proj);
y = H.apply(vol_coeff);
sigma = sqrt(var(y(:))/10);
y = y + sigma_noise(sig_noise_index) * randn(size(y));
size_im = [size(y,1), size(y,2)];

% define operators
define_operators;

% intializing the volume
% one can replace vol_init with any other 3D ab-initio model
vol_init = imgaussfilt3(vol_gt, 4);
[vol_coeff_init, vol_init_est, kbwf_recon_kernel] = expand_vol_KBWF(vol_init, kbwf_recon);
plot_vol_results(vol_init, vol_init_est, 'Init volume', 'Expanded init volume')

% initializing angles and shifts
angles_init = angles + sigma_angle*(rand(size(angles))-0.5)*2;
shifts_init = zeros(size(shifts));

known_angle = true;
approx_angle = true;
joint_angle_vol = true;

%% oracle baseline: refine based on ground truth (GT) angles
if known_angle
    H_true = LinOpPBTShift(size(vol_coeff), angles, shifts, 0, kbwf_recon);
    [vol_coeff_true, snr_true_evol, iter_true_evol, vol_true_recon_iter] = ...
        ADMM_solver(10, LS * H_true, Fn, Hn, [rho_n_final, 1e4], zeros(size(vol_coeff)), vol_coeff_init, 1); % vol_init
    vol_true = convn(vol_coeff_true, kbwf_recon_kernel, 'same');
    plot_vol_results(vol_gt, vol_true, 'GT volume', 'Recon. volume')
end

%% approx baseline: refine based on erroneous angles
if approx_angle
    H_approx = LinOpPBTShift(size(vol_coeff), angles_init, shifts_init, 0, kbwf_recon);
    [vol_coeff_approx, snr_approx_evol, iter_approx_evol, vol_approx_recon_iter] = ...
        ADMM_solver(10, LS * H_approx, Fn, Hn, [rho_n_final, 1e4], zeros(size(vol_coeff)), vol_coeff_init, 1);
    vol_approx = convn(vol_coeff_approx, kbwf_recon_kernel, 'same');
    plot_vol_results(vol_gt, vol_approx, 'GT volume', 'Recon. volume')
end

%% joint optimization of the volume and the angles
if joint_angle_vol
    vol_mask = ones(size(vol_gt));
    struct_params;
    [ vol_rec_final, angles_rec_iter, shifts_rec_iter, final_SNR_iter, vol_recon_cell, ~, ~] ...
        = alternating_angle_vol(params);
    vol_joint = convn(vol_rec_final, kbwf_recon_kernel, 'same');
    plot_vol_results(vol_gt, vol_joint, 'GT volume', 'Recon. volume')
end

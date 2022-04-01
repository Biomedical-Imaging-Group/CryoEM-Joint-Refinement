clear all
clc

rng(1);

%% loading the required parameters
param_exp;
num_samples = 5000; % # of projs

%% optimization setting and definition
optim_config;

% generating the volume
vol_name = 'emd_3400';
vol_gt = prepare_volume(vol_name, 1);
vol_sz = size(vol_gt, 1);

[vol_coeff, vol_est, kbwf_kernel] = expand_vol_KBWF(vol_gt, kbwf_proj);
% vol_coeff = vol_gt;
plot_vol_results(vol_gt, vol_est, 'GT volume', 'Expanded volume')

% generate the angles and the shifts
sigma_shiftX = 4;
sigma_shiftY = 4;
[rot, tilt, psi] = generateEquidistribRandomProjAngles(num_samples);
shifts = ShiftGen(sigma_shiftX, sigma_shiftY, num_samples); % used to create the shifts
angles = [rot',tilt',psi'];

%% Range of the parameters of the experiment
% different noise levels
sigma_noise = 10.^(linspace(-1, 2, 5));
sig_noise_index = 1; %4
sigma_angle = 0.7;

%% generating the forward model and other operators
H = LinOpPBTShift(size(vol_coeff), angles, shifts, 1, kbwf_proj);
y_clean = H.apply(vol_coeff);
y = y_clean + sigma_noise(sig_noise_index) * randn(size(y_clean));
size_im = [size(y,1),size(y,2)];

% define operators
define_operators;

% intializing the volume
vol_init = imgaussfilt3(vol_gt, 4); %4 for 46
% vol_init = vol_init / max(vol_init(:));
% vol_init = lp_filter_volume(vol_gt, 0.06);
[vol_coeff_init, vol_init_est, kbwf_recon_kernel] = expand_vol_KBWF(vol_init, kbwf_recon);
% vol_coeff_init = vol_init;
plot_vol_results(vol_init, vol_init_est, 'Init volume', 'Expanded init volume')

% initializing angles
angles_init = angles + sigma_angle*(rand(size(angles))-0.5)*2;

% initializing shifts
shifts_init = zeros(size(shifts));

known_angle = true;
approx_angle = true;
joint_angle_vol = true;

%% oracle baseline: refine based on ground truth (GT) angles
if known_angle
    H_true = LinOpPBTShift(size(vol_coeff), angles, shifts, 0, kbwf_recon);
    [vol_coeff_true, snr_true_evol, iter_true_evol, vol_true_recon_iter] = ...
        ADMM_solver(10, LS*H_true, Fn, Hn, [rho_n_final,10000], zeros(size(vol_coeff)), vol_coeff_init, 1); % vol_init
    vol_true = convn(vol_coeff_true, kbwf_recon_kernel, 'same');
    plot_vol_results(vol_gt, vol_true, 'GT volume', 'Recon. volume')
end

%% approx baseline: refine based on erroneous angles
if approx_angle
    H_approx = LinOpPBTShift(size(vol_coeff), angles_init, shifts_init, 0, kbwf_recon);
    [vol_coeff_approx, snr_approx_evol, iter_approx_evol, vol_approx_recon_iter] = ...
        ADMM_solver(20, LS*H_approx, Fn, Hn, [rho_n_final,10000], vol_coeff, vol_coeff_init, 1);
    vol_approx = convn(vol_coeff_approx, kbwf_recon_kernel, 'same');
    plot_vol_results(vol_gt, vol_approx, 'GT volume', 'Recon. volume')
end

%% joint optimization of the volume and the angles
if joint_angle_vol
    vol_mask = ones(size(vol_gt));
    struct_params;
    [ vol_rec_final, angles_rec_iter, shifts_rec_iter, final_SNR_iter, vol_recon_cell, ~, ~] ...
        = alternating_angle_vol2(params);
end

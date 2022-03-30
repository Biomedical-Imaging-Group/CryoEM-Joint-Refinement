% define the variables input to the joint reconstruction pipeline
% Mona Zehni, 2022

params = struct;
params.max_iter = max_iter;
params.vol_init = vol_coeff_init;
params.angles_init = angles_init;
params.shifts_init = shifts_init;
params.projs = y;

params.LS = LS;
params.R1 = R1;
params.rho_n_init = rho_n_init;
params.rho_n_final = rho_n_final;
params.max_iter_ADMM = max_iter_ADMM;
params.grad_num = grad_num;
params.grad_step_rot = grad_step_rot;
params.grad_step_tilt = grad_step_tilt;
params.grad_step_shift = grad_step_shift;
params.alpha = alpha;
params.beta = beta;
params.num_scale_stages = num_scale_stages;
params.lamb_init = lamb1_init;
params.lamb_final = lamb1_final;
params.vol_mask = vol_mask;
params.DL_proj = 1;

% only used to compare the estimated values against
params.angles = angles;
params.shifts = shifts;

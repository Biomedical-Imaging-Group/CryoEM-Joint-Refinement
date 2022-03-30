%the set of parameters
% Mona Zehni, 2018

scale1_rho_init = [2e3, 2e3, 2e3, 7e4, 3e5];
scale1_lamb_init = [100, 100, 100, 1e3, 4e3];

scale1_rho_final = [2e2, 2e2, 2e2, 3.2e4, 3e5];
scale1_lamb_final = [10, 10, 1e1, 500, 4e3];

kbwf_proj.a = 2;
kbwf_proj.m = 2;
kbwf_proj.alpha = 19;
kbwf_proj.scale = 1;
kbwf_proj.Ti = 1;
kbwf_proj.Tp = 1;

kbwf_recon.a = 4;
kbwf_recon.m = 2;
kbwf_recon.alpha = 19;
kbwf_recon.scale = 1;
kbwf_recon.Ti = 1;
kbwf_recon.Tp = 1;
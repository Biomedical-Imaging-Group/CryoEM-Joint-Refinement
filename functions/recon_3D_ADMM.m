function [vol, final_SNR, final_fval, final_iternum] = ...
    recon_3D_ADMM(LS, H, Fn, Hn, rho_n, vol_true, vol_init, max_iter)
%Performs several ADMM steps
%Input:
%       LS, H, Fn, Hn, rho_n: the operators and the parameters of
%       optimization problem
%       vol_true: the true volume, only used for computing the gradientds
%       at each step
%       vol_init: the initial volume to start the optimizations with
%       max_iter: maximum number of iterations (or the number of ADMM iterations)
%Output:
%       vol: the updated volume after taking several steps of ADMM
%       final_SNR: the SNR evolution throughout the ADMM iterations
%       final_fval:the cost value evolution throughout the ADMM iterations
%       final_iternum: the evolution of the iteration number throughout the
%       ADMM steps
% Mona Zehni, July 2018

% eps_abs = 1.5e-2;
% eps_rel = 2e-3;

% eps_abs = 1e-1;
% eps_rel = 2e-3;

ADMM=OptiADMM(LS*H,Fn,Hn,rho_n);
% ADMM.OutOp=OutputOptiSNR(0,vol_true,1);

ADMM.ItUpOut=1;
ADMM.maxiter=max_iter; 
ADMM.run(vol_init);  
vol = ADMM.xopt;

final_SNR = 0; %ADMM.OutOp.evolsnr;
final_fval = ADMM.OutOp.evolcost;
final_iternum = ADMM.OutOp.iternum;

end
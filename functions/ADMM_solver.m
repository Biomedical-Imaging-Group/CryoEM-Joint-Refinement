function [final_vol, snr_evol, iter_evol, vol_recon_iter] = ADMM_solver(maxIter, F0, Fn, Hn, rho_n, vol_true, vol_init, scaleRatio)
%Solves and defines the regularized optimization problem using ADMM
%Input:
%       maxIter: Maximum number of used iterations
%       F0, Fn, Hn, rho_n: the operators used in the definition of the
%       inverse problem
%       vol_true: The true volume, used in order to compute SNR in
%       different stages
%       vol_init: the initial volume to start the ADMM process with
%       scaleRatio: the scaling of the recovery, 1 is the fines scale
%Output:
%       final_vol: the final recovered volume after using ADMM
%       snr_evol: the SNR evolution in the ADMM process
%       iter_vol: the corresponding iteration number to the saved snr_evol
% Mona Zehni, July 2018

eps_abs = 1e-1;
eps_rel = 2e-3;

ADMM=OptiADMM(F0,Fn,Hn,rho_n);

% the SNR should be computed wrt a downsampled version of the volume if the
% scaling ratio is not one
% vol_dl = downSample(vol_true, scaleRatio);
% ADMM.OutOp=OutputOpti(0,vol_dl,1);
% ADMM.CvOp=TestCvgADMM(eps_abs,eps_rel);

ADMM.ItUpOut=1;             
ADMM.maxiter=maxIter; 
ADMM.run(vol_init);

final_vol = ADMM.OutOp.evolxopt{end};
snr_evol = ADMM.OutOp.evolsnr;
iter_evol = ADMM.OutOp.iternum;

vol_recon_iter = ADMM.OutOp.evolxopt;
end


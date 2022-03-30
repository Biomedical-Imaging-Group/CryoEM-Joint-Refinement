% define the operators
G = LinOpGrad(size(vol_coeff)); % gradient operator
R1 = CostL1(G.sizeout, zeros(G.sizeout)); % total variation regularizer
Hn = {G, LinOpIdentity(size(vol_coeff))};
R_pos = CostNonNeg(size(vol_coeff)); % positivity constraint

LS = CostL2([],y); % least squares formulation
% starting and ending rhos and lambdas used by the alternating model
rho_n_init = scale1_rho_init(sig_noise_index)/1;
lamb1_init  = scale1_lamb_init(sig_noise_index)/50; %/1

rho_n_final = scale1_rho_final(sig_noise_index)/50;
lamb1_final  = scale1_lamb_final(sig_noise_index)/500;

%% oracle baseline: true angles are known
Fn = {lamb1_final*R1, R_pos};
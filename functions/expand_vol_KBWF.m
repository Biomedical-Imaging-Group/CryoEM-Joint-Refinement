function [c, vol_est, kbwf] = expand_vol_KBWF(vol, kbwf_config)
%% Expands volume on KBWF basis
%input vol: the volume to be expanded
%input kbwf_config: parameters of the KBWF

%output c: expansion coefficients
%output vol_est: estimated volume based on the expansion
%output kbwf: kbwf kernel
% Mona Zehni, 2022

supp = kbwf_config.a;
m = kbwf_config.m;
alpha = kbwf_config.alpha;
scale = kbwf_config.scale;

vol_sz = size(vol, 1); % assume the volume is cubic

n = floor(vol_sz/2);
if mod(vol_sz, 2) == 0
    grid = [-n:n - 1];
else
    grid = [-n:n];
end

[grid_x, grid_y, grid_z] = meshgrid(grid, grid, grid);
s = sqrt(grid_x.^2 + grid_y.^2 + grid_z.^2);
tmp = sqrt(max(1 - (s/supp).^2, 0));
kbwf = 1 ./ besseli(m, alpha) .* tmp.^(m) .* besseli(m, alpha*tmp);

% form an optimization problem to find the expansion coeffs
H = LinOpConv(fftn(kbwf));
L2 = CostL2(H.sizeout, vol);  
F0 = L2*H;
G = LinOpGrad(size(vol)); % gradient operator
R1 = CostL1(G.sizeout, zeros(G.sizeout)); % total variation regularizer
Hn = {G, LinOpIdentity(size(vol))};
R_pos = CostNonNeg(size(vol)); % positivity constraint
Fn = {0.001*R1, R_pos};
ADMM=OptiADMM(F0,Fn,Hn,[1, 1]);

ADMM.ItUpOut=50;             
ADMM.maxiter=200; 
ADMM.run(zeros(size(vol)));

c = ADMM.xopt;
c = fftshift(c);

vol_est = convn(c, kbwf, 'same');
end
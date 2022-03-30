% optimization config
% Mona Zehni, 2019

grad_num = 3; 
max_iter_ADMM = 5; 
maxIter = 200; % for the other optimizations
grad_step_rot = 10^-5;%1e-7; 
grad_step_tilt = 10^-5;%1e-7; 
grad_step_shift = 1e-4; %1e-5

% parameters of the line search
alpha = 0; 
beta = 0.25; 

max_iter = 10; % total number of iterations for the alternating approach
num_scale_stages = 1; % number of scale stages in the refinement
% main script for implementing the GPLFM on an MDOF cantilever model with acceleration-only measurements. 
% author: Joanna Zou
% version: 0.1.0, @2023

clear; close all; clc;
addpath(pwd, genpath(pwd));


%% FOR EDITING: user-defined parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A) mechanical model 
ndof = 10;                       % number of stories (degrees of freedom)
mass = 200;                      % mass of each floor (kg)
k = 5*10^3;                      % spring constant (N/m)   

% B) measured response 
Sa_dof = 1:2:10;                 % measured dofs of accelerations

% C) simulation and filtering
dt_f = 0.001;                    % time step for simulating force 
fs_f = 1/dt_f;                   % simulation frequency
dt_n = 0.05;                     % time step for downsampling simulated response; lower for better accuracy
dt = 0.01;                       % time step for filtering; lower for better accuracy
fs = 1/dt;                       % filtering frequency
T = 1000;                        % final time
% See file 'define_load_synthetic.m' to change input parameters

% D) state-space models
nmodes = 3;                      % number of modes for model order reduction
sigQ = 10^(-16);                 % process noise variance
sigR = 10^(-4);                  % measurement noise variance

% E) HP optimization 
niter = 1;                       % number of random seeds
hp_bd = [0.01, 100;              % alpha lower bd. and upper bd.
         0.0001, 10];            % ls lower bd. and upper bd.



%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A) create model
[M,K,C] = model_ShearBldg(ndof, mass, k);


%% B) define measurements
Sa = zeros(length(Sa_dof),ndof);
for d = 1:length(Sa_dof)
    Sa(d,Sa_dof(d)) = 1;
end


%% C) simulate ground truth response

% define load
define_load_synthetic;

% simulate response by average Newmark acceleration method
x_0 = zeros(2*ndof,1);
Mdl_sim = ss_MDOF(M,K,C,Sa,Sp,dt_f);
[t, x_sim, y_sim, ~, ~] = rs_Newmark(Mdl_sim, loads, 1/2, 1/4, x_0, sigQ, sigR);
acc_n = Sa*y_sim;

% downsample the simulated response (for comparison)
[~, x_sim_n] = downsample(t, x_sim, dt);
[~, y_sim_n] = downsample(t, y_sim, dt);

% downsample the measurement signal
[t_n, acc_n] = downsample(t, acc_n, dt_n);


%% D) create state-space models

% modally reduced-order model for filtering
Mdl = ss_MDOF_modal(nmodes,M,K,C,Sa,dt);
Qx = eye(Mdl.nx)*sigQ;       % process noise for x       
R = eye(Mdl.ny)*sigR;        % measurement noise 


%% E) GP hyperparameter inference 
time = 0:dt:T;                   % time vector for filtering
nsteps = length(time);           % number of time steps

% reduced order model
opt = hp_Optimizer(niter, Mdl, hp_bd, acc_n);
hpopt = opt.hp_opt(:,opt.idx);


%% F) sequential inference of response states
                                        
% reduced order model
[z_r, cov_r, a_r, cova_r, MdlA, ~] = id_GPLFM(t_n, acc_n, time, Mdl, @kernel_Matern52, hpopt, opt.x_0, R, 0);


% convert modal states to original states: 

% displacement/velocity states
x_r = Mdl.Transf*z_r(1:Mdl.nx,:);

% recover modal forces     
d_r = Mdl.nx + 1;
f_r = z_r(d_r:3:end,:);         

% transform reference force into modal coordinates
f_true = Mdl.Phi_r' * Sp * loads.f;
[~, f_true] = downsample(t, f_true, dt);


%% make plots
plot_GPLFM_Modal;


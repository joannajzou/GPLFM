function [knl, K] = kernel_Matern52(hp, dt, xa, xb)
% returns covariance matrix for a GP with Matern 5/2 kernel

% INPUT:
% hp (1d vector) = list of hyperparameters for the kernel.
%   hp(1) = alpha (vertical length scale)
%   hp(2) = ls (horizontal length scale)
% dt = time step for discretized matrices
% xa (na x 1) [OPT] = list of points to condition on ("observations")
% xb (nb x 1) [OPT] = list of linearly spaced points to predict ("tests")

% OUTPUT: 
% kernel (struct) = contains state-space matrices of the kernel
% K (na x nb) = covariance matrix defined by the kernel



    knl = struct;

    % kernel hyperparameters
    nu = 5/2;
    alpha = hp(1);
    ls = hp(2);
    knl.hp = hp;
    
    lambda = sqrt(5)/ls;
    knl.q = (400*sqrt(5)*alpha^2) / (3*ls^5); % spectral density

    
    % continuous state-space matrices
    knl.F_c = zeros(3,3);
    knl.F_c(1:2,2:3) = eye(2);
    knl.F_c(3,1) = -lambda^3;
    knl.F_c(3,2) = -3*lambda^2;
    knl.F_c(3,3) = -3*lambda;

    knl.L_c = [0; 0; 1];
    knl.H_c = [1, 0, 0];
    
    % discrete state-space matrices
    knl.nz = size(knl.F_c,1);
    knl.F = expm(knl.F_c*dt);
    knl.L = (knl.F - eye(knl.nz))*inv(knl.F_c)*knl.L_c;
    knl.H = knl.H_c;

    % solve for steady-state covariance from the continuous Lyapunov eqn
    knl.P = lyap(knl.F_c, knl.L_c*knl.q*knl.L_c');
    % set up prior for the state z
    knl.Q = knl.P - knl.F*knl.P*knl.F.';
    
    knl.sigma_w = sqrt(knl.q/dt);
    
    
    
    if ~exist('xa','var') & ~exist('xb','var')
        K = [];
    else
        % if xa, xb supplied, then find covariance matrix
        K = zeros(length(xa), length(xb));

        for i = 1:length(xa)
            for j = 1:length(xb)
                tau = sqrt((xa(i) - xb(j)).^2);
                K(i,j) = (alpha^2).*(1 + (sqrt(5).*tau./ls) + (5.*tau.^2./(3*ls^2)))*exp(-sqrt(5).*tau./ls);
            end
        end
    end
    
end    
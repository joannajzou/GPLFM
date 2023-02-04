function [t, x_sim, y_sim, Q, R] = rs_Newmark(mdl, loads, gamma, beta, x_0, sigQ, sigR)
% uses Newmark method to simulate "ground truth" response time history for a linear system
% average acceleration method: gamma=1/2, beta=1/4 (stable for any time step)
% linear acceleration method: gamma=1/2, beta=1/6 (stable for dt/T1<=0.0551)

% INPUT:
% mdl (struct) = contains model and measurement matrices
% loads (struct) = contains force time series and properties
% gamma (float) = gamma parameter for Newmark method (typ. set to 1/2)
% beta (float) = beta parameter for Newmark method (1/4 or 1/6)
% x_0 (nx-by-1) = initial condition for states x
% sigQ (float) = variance of process noise on states x
% sigR (float) = variance of measurement noise on states y

% OUTPUT:
% t (1-by-nsteps) = time vector
% x_sim (nx-by-nsteps) = state time series at each dof
% y_sim (ny-by-nsteps) = response time series at each dof
% Q (nx-by-nx) = covariance matrix of process noise
% R (ndof-by-ndof) = covariance matrix of measurement noise


    nsteps = size(loads.f,2);
    dt = mdl.dt;
    t = 0:dt:(nsteps-1)*dt;
    
    
    % initialize vectors and set initial conditions
    dis = zeros(mdl.ndof, nsteps);
    dis(:,1) = x_0(1:mdl.ndof);
    vel = zeros(mdl.ndof, nsteps);
    vel(:,1) = x_0(mdl.ndof+1:end);
    acc = zeros(mdl.ndof, nsteps);
    
    
    % initial calculations
    acc(:,1) = inv(mdl.M)*(mdl.S_p*loads.f(:,1) - mdl.C*vel(:,1) - mdl.K*dis(:,1));
    a1 = (1/(beta*dt^2))*mdl.M + (gamma/(beta*dt))*mdl.C;
    a2 = (1/(beta*dt))*mdl.M + ((gamma/beta)-1)*mdl.C;
    a3 = (1/(2*beta)-1)*mdl.M + dt*((gamma/(2*beta))-1)*mdl.C;
    Khat = mdl.K + a1;
    

    % simulation
    for k = 1:nsteps-1
        Phat = mdl.S_p*loads.f(:,k+1) + a1*dis(:,k) + a2*vel(:,k) + a3*acc(:,k);
        dis(:,k+1) = inv(Khat)*Phat;
        vel(:,k+1) = (gamma/(beta*dt))*(dis(:,k+1) - dis(:,k)) + (1-(gamma/beta))*vel(:,k) + dt*(1-(gamma/(2*beta)))*acc(:,k);
        acc(:,k+1) = (1/(beta*dt^2))*(dis(:,k+1) - dis(:,k)) - (1/(beta*dt))*vel(:,k) - ((1/(2*beta))-1)*acc(:,k);
    end
    

    x_sim0 = [dis; vel];
    y_sim0 = acc;


    % simulate process noise
    Q = sigQ * eye(mdl.nx);
    w = mvnrnd(zeros(mdl.nx,1), Q, nsteps)';
    x_sim = x_sim0 + w;
    

    % simulate measurement noise
    R = sigR * eye(mdl.nye); 
    v = mvnrnd(zeros(mdl.nye,1), R, nsteps)';
    y_sim = y_sim0 + v;
    

end





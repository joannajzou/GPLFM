function [mu_ab, cov_ab, out, cov_out, mdlA, R] = id_GPLFM(x1, y1, x2, mdl, kernel_func, hp, mu_0, R_0, tuneflag) 
% performs joint input-state estimation with the GPLFM via Kalman filtering

% INPUT:
% x1 (1-by-nt) = input (time) vector corresponding to measurements
% y1 (ny-by-nt) = measurements vector 
% x2 (1-by-nt_n) = input (time) vector corresponding to filtering predictions
% mdl (struct) = contains model and measurement SS matrices
% kernel_func (func) = function generating kernel matrices
% hp (nhp-by-1) = GP hyperparameters
% mu_0 (nx-by-1) = prior mean on the states
% R_0 (ny-by-ny) = initial value of measurement noise covariance
% tuneflag(1 or 0) = if 0, then measurement noise tuning is not performed


% OUTPUT:
% mu_ab (nx-by-nsteps) = posterior mean state estimate
% cov_ab (1-by-nsteps cell) = cell array which contains the nx-by-nx posterior covariance matrices on the state
% out (nx-by-nsteps) = posterior mean estimate on outputs
% cov_out (1-by-nsteps cell) = cell array which contains the nx-by-nx posterior covariance matrices on the outputs
% mdlA (struct) = augmented SS matrices
% R (ny-by-ny) = final value of measurement noise covariance after tuning


    % allocate vectors
    nsteps = size(x2,2);
    dt_n = x1(2)-x1(1);
    dt_f = x2(2)-x2(1);
                         
    
    % set kernel for each force
    knlstruct = {};
    for p = 1:mdl.np
        [knlstruct{p}, ~] = kernel_func(hp, dt_f);
    end
 
    % create augmented state-space model
    mdlA = ss_Augmented(mdl, knlstruct, dt_f);
                                                         

    % initial conditions
    mu_af = zeros(mdlA.na, nsteps);                      % forward pass mean
    mu_af(1:mdl.nx,1) = mu_0;                            % assume zero initial conditions for force
    cov_af = {};                                         % forward pass covariance

    out = zeros(mdl.nye, nsteps);                        % output states (acc.)
    cov_out = {};                                        % covariance on output states
    
    cov_out_0 = mdlA.H_ae*mdlA.P_0*mdlA.H_ae';
    cov_y_0 = mdlA.H_ad*mdlA.P_0*mdlA.H_ad';             % for computing the prior

    for n = 1:nsteps
        cov_af{n} = mdlA.P_0;  
        cov_out{n} = cov_out_0;
    end

    mu_ab = mu_af;                                       % backward pass mean
    cov_ab = cov_af;                                     % backward pass covariance
    
    mu_pred = mu_af;                                     % prediction mean
    cov_pred = cov_af;                                   % prediction covariance


    % measurement noise tuning parameters
    e_k = zeros(mdl.ny, nsteps);                         % residual
    R = R_0;                                             % initialize noise covariance
    r_mag = R(1,1);
    tol1 = 0.1;
    tol2 = 1e-9;
    conv = 0;
    numiter = 0;
    

    % find indices of measurements
    ind = ismember(round(x2,4),round(x1,4));
    ind = find(ind==1);
    y2 = zeros(size(y1,1),nsteps);
    y2(:,ind) = y1;                                      % extend y1 vector to 1-by-nsteps


    % evaluate fit of HP
    y1_emp = y1(find(y1(:,1)~=0),:);
    cov_y_emp = cov(y1_emp');
    plot_EvalGPPrior(x1, y1, cov_y_emp, cov_y_0);
    

    % start loop
    k = 1;
    count = 1;
    
    %% forward pass: filtering
    
   while conv == 0
       for k = 2:nsteps
            % Predict step
            mu_pred(:,k) = mdlA.F_ad*mu_af(:,k-1);
            cov_pred{k} = mdlA.F_ad*cov_af{k-1}*mdlA.F_ad' + mdlA.Q_ad;
            % Update step
            if any(ismember(k, ind))
                e_k(:,k) = y2(:,k) - mdlA.H_ad*mu_pred(:,k);
                S_k = mdlA.H_ad*cov_pred{k}*mdlA.H_ad' + R;
                K_k = cov_pred{k}*mdlA.H_ad'*inv(S_k);
                mu_af(:,k) = mu_pred(:,k) + K_k*e_k(:,k);
                cov_af{k} = (eye(mdlA.na) - K_k*mdlA.H_ad)*cov_pred{k}*(eye(mdlA.na) - K_k*mdlA.H_ad)' + K_k*R*K_k';


            else
                % if no observation is made, no update
                mu_af(:,k) = mu_pred(:,k);
                cov_af{k} = cov_pred{k};

            end
            out(:,k) = mdlA.H_ae*mu_af(:,k);
            cov_out{k} = mdlA.H_ae*cov_af{k}*mdlA.H_ae';
            
       end
       conv = 1; 

       if tuneflag == 1
            % check noise covariance
            r = var(e_k');
            R_f = median(nonzeros(r))*diag(r==0) + diag(r);
            R_error = diag(abs(R_f - R))./max(diag(R));
            if max(R_error) < tol1 % & max(r) < tol2
                conv = 1;
            else
                R = R_f;
            end
            numiter = numiter + 1
       else
           conv = 1;
   end
        
    
    %% backward pass: smoothing

    % initialize
    mu_ab(:,nsteps) = mu_af(:,nsteps);
    cov_ab{nsteps} = cov_af{nsteps};
    

    for k = flip(1:nsteps-1)% iterate starting from last step
        % smoothing update
        Ks = cov_af{k}*mdlA.F_ad'*inv(cov_pred{k+1});
        [a, MSGID] = lastwarn();
        warning('off', MSGID);
        mu_ab(:,k) = mu_af(:,k) + Ks*(mu_ab(:,k+1)-mu_pred(:,k+1));
        cov_ab{k} = cov_af{k} + Ks*(cov_ab{k+1} - cov_pred{k+1})*Ks';

        out(:,k) = mdlA.H_ae*mu_ab(:,k);
        cov_out{k} = mdlA.H_ae*cov_ab{k}*mdlA.H_ae';
    end
    

%     % or else skip smoothing
%     mu_ab = mu_af;
%     cov_ab = cov_af;

    
end

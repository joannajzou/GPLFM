function dist = hp_MinDist(hp_norm, y, mdl, opt, kernel_func)
% objective function for fitting GP hyperparameters by minimizing variance 

% INPUT:
% hp_norm (nhp x 1) = vector of GP hyperparameters, normalized to range from 0 to 1
% y (ny-by-nsteps) = measurement time series
% mdl (struct) = contains model SS matrices
% opt (struct) = contains optimization settings
% kernel_func (func) = name of function for computing the kernel SS matrices

% OUTPUT:
% dist (scalar) = Hellinger dist. between modeled and empirical distributions


    % rescale to original coordinates
    hp = hp_rescale(hp_norm, [], opt.bounds); 
    
    % define kernel for each force
    knlstruct = {};
    for p = 1:mdl.np
        [knlstruct{p}, ~] = kernel_func(hp, mdl.dt);
    end

    % construct augmented state-space model
    mdlA = ss_Augmented(mdl, knlstruct, mdl.dt);
                                                        
    % compute modeled covariance matrix on measured states
    cov_y_hp = mdlA.H_ad*mdlA.P_0*mdlA.H_ad';
    
    % compute empirical covariance matrix on measured states
    y_emp = y(find(y(:,1)~=0),:);
    ny = size(y_emp,1);
    cov_y_emp = cov(y_emp');
    
    % compute Hellinger distance between two distributions 
    Stat_emp = [zeros(ny,1),sqrt(abs(diag(cov_y_emp)))];    % empirical dist.
    Stat_hp = [zeros(ny,1),sqrt(abs(diag(cov_y_hp)))];      % modeled dist.
        
    dist = val_Hellinger(Stat_emp,Stat_hp,500);
    
end

function opt = hp_Optimizer(niter, mdl, bds, data)
% solve for GP hyperparameters

% INPUT:
% niter (int) = number of random initializations
% mdl (struct) = contains model and measurement matrices
% bds (nhp-by-2) = lower and upper bounds of search domain
% data (ny-by-nsteps) = measurement data 

% OUTPUT:
% opt (struct) = contains optimization constraints and results


    opt = struct;                           % initialize opt struct
    opt.nhp = 2;                            % number of HPs
    opt.x_0 = zeros(mdl.nx, 1);             % prior mean
    opt.bounds = bds;                       % hp ll and ul
    opt.hp0 = zeros(opt.nhp,niter);         % initial guess for HPs
    opt.hp0_norm = zeros(opt.nhp,niter);    % initial guess for HPs
    opt.hp_opt = zeros(opt.nhp,niter);      % optimal HPs
    opt.fval = zeros(1,niter);              % func. eval at optimum
    opt.hp_hist = {};                       % iterates of HPs
    opt.fval_hist = {};                     % iterates of func. eval
    
    % objective function
    opt.objective = @(theta) hp_MinDist(theta, data, mdl, opt, @kernel_Matern52);
    
    
    for iter = 1:niter
        % randomly sample initial point
        hp0_norm = rand(opt.nhp,1);         % normalized on scale of 0 to 1
        opt.hp0_norm(:,iter) = hp0_norm;
        opt.hp0(:,iter) = hp_rescale(hp0_norm, [], opt.bounds); % in true coordinates
        
        % run optimization
        [hp_opt_norm, fval, opt.hp_hist{iter}] = hp_fmincon(hp0_norm, opt.objective,...
                                              [],[],[],[],[0;0],[1;1]);           
        opt.hp_opt(:,iter) = hp_rescale(hp_opt_norm, [], opt.bounds);
        opt.fval(iter) = fval;
        
        % find func. eval over iterates
        opt.fval_hist{iter} = [];
        for i = 1:size(opt.hp_hist{iter},2)
            % hpi = hp_rescale(opt.hp_hist{iter}(:,i), [], opt.bounds);
            opt.fval_hist{iter} = [opt.fval_hist{iter}, opt.objective(opt.hp_hist{iter}(:,i))];
        end
    end
    
    
    % find minimum func eval across all iterations
    [opt.fmin, opt.idx] = min(opt.fval);   


end



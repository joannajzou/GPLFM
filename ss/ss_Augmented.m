function mdlA = ss_Augmented(mdl, knl, dt)
% creates augmented state-space matrices of the GPLFM

% INPUT:
% mdl (struct) = contains model and measurement SS matrices
% knl (np cell) = contains structs with kernel SS matrices corresponding
%                 to each of np inputs
% dt (float) = discretization time step

% OUTPUT:
% mdlA (struct) = contains augmented model and measurement SS matrices 


    mdlA = struct;
    
    % compute augmented dimension (function of number of inputs)
    dim = 0;
    for m = 1:mdl.np
        dim = dim + knl{m}.nz;
    end

    % create continuous augmented matrices
    mdlA.F_cs = zeros(dim,dim);
    mdlA.B_cs = zeros(mdl.nx,dim);
    mdlA.J_cs = zeros(mdl.ny,dim);
    mdlA.J_cse = zeros(mdl.nye,dim);
    mdlA.Q_c = zeros(dim,dim);
    j = 1;
    for m = 1:mdl.nm
        block = knl{m}.nz-1;
        mdlA.F_cs(j:j+block,j:j+block) = knl{m}.F_c;
        mdlA.B_cs(:,j:j+block) = mdl.B_c(:,m)*knl{m}.H_c;
        mdlA.J_cs(:,j:j+block) = mdl.J_c(:,m)*knl{m}.H_c;
        mdlA.J_cse(:,j:j+block) = mdl.J_e(:,m)*knl{m}.H_c;
        mdlA.Q_c(j:j+block,j:j+block) = knl{m}.L_c*knl{m}.q*knl{m}.L_c';
        j = j+knl{m}.nz;
    end

    mdlA.F_ac = [mdl.A_c mdlA.B_cs; ...
                zeros(dim, mdl.nx) mdlA.F_cs];
    mdlA.H_ac = [mdl.G_c mdlA.J_cs];
    mdlA.H_ae = [mdl.G_e mdlA.J_cse];

    mdlA.na = size(mdlA.F_ac,1);


    % discretized augmented matrices
    mdlA.F_ad = expm(mdlA.F_ac*dt);
    mdlA.H_ad = mdlA.H_ac;

    
    % process noise covariance matrix - solve with Lyapunov equation
    Qx = zeros(mdl.nx); 
    mdlA.Q_ac = [Qx, zeros(mdl.nx, dim);... 
                zeros(dim, mdl.nx), mdlA.Q_c];
    
    mdlA.P_0 = lyap(mdlA.F_ac, mdlA.Q_ac);
    
    mdlA.Q_ad = mdlA.P_0 - mdlA.F_ad*mdlA.P_0*mdlA.F_ad';
    

end 
    
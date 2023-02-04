function mdl = ss_MDOF_modal(n,M,K,C,S_acc,dt,maxvals)
% returns modally reduced state-space matrices describing the model, measurements, and input

% INPUT:
% n (int) = number of modes to include in the reduced-order model
% M (ndof-by-ndof) = model mass matrix 
% K (ndof-by-ndof) = model stiffness matrix 
% C (ndof-by-ndof) = model damping matrix 
% S_acc (ndof-by-ndof) = indicator matrix of measured accelerations
% dt (float) = time step
% maxvals (vector) [OPT] = vector of maximum values, for normalizing outputs

% OUTPUT:
% mdl (struct) = model with the system matrices as added properties
%   A_c (nx-by-nx) = continuous state matrix 
%   B_c (nx-by-np) = continuous measurement matrix 
%   G_c (ny-by-ny) = continuous output influence matrix 
%   J_c (ny-by-np) = continuous direct transmission matrix
%   A (nx-by-nx) = discrete state matrix 
%   B (nx-by-np) = discrete measurement matrix 
%   G (ny-by-ny) = discrete output influence matrix 
%   J (ny-by-np) = discrete direct transmission matrix


    mdl = struct;   
    
    % save inputs in mdl struct
    mdl.M = M;                                     
    mdl.K = K;                                      
    mdl.C = C;   
    mdl.Sa = S_acc;
    mdl.Sp = eye(n); % for modal force estimation
    mdl.dt = dt;
    
    % compute modal properties
    [Phi, Omega2] = eig(K,M);
    w = sqrt(diag(Omega2));
    f = w/(2*pi);
    Gamma = diag(diag(Phi'*C*Phi)); % double "diag" ensures off-diagonal terms are zero
    zeta = diag(Gamma)/2./w;
    
    % compute reduced-order matrices
    mdl.Phi_r = Phi(:,1:n);
    mdl.Omega2_r = Omega2(1:n,1:n);
    mdl.Gamma_r = Gamma(1:n,1:n);
    mdl.Transf = [mdl.Phi_r zeros(size(mdl.Phi_r));...
                  zeros(size(mdl.Phi_r)) mdl.Phi_r];
              mdl.f = f(1:n)';
    mdl.zeta = zeta(1:n);
              
    % define reduced-order model dimensions
    mdl.ndof_orig = size(mdl.M,1);
    mdl.ndof = n;                                   % number of dofs
    mdl.nm = n;                                     % number of modes
    mdl.nx = 2*n;                                   % number of states
    mdl.ny = size(S_acc,1);                         % number of outputs
    mdl.nye = mdl.ndof_orig;                        % number of outputs to predict
    mdl.np = n;                                     % number of external (modal) forces
    
    % modal force 
    mdl.P = [eye(mdl.np)];
    
    % compute continuous-time system matrices
    mdl.A_c = [zeros(mdl.ndof) eye(mdl.ndof);...
              -mdl.Omega2_r -mdl.Gamma_r];
      
    mdl.B_c = [zeros(mdl.ndof,mdl.np); mdl.P];
       
    mdl.G_c = [-S_acc*mdl.Phi_r*mdl.Omega2_r -S_acc*mdl.Phi_r*mdl.Gamma_r];
    mdl.J_c = [S_acc*mdl.Phi_r*mdl.P];
    
    % for estimating all acceleration outputs
    mdl.G_e = [-mdl.Phi_r*mdl.Omega2_r -mdl.Phi_r*mdl.Gamma_r];
    mdl.J_e = [mdl.Phi_r*mdl.P];
    
    
    % if specified, scale G_c to normalize outputs 
    if exist('maxvals','var')
        mdl.G_c = ((mdl.G_c')./(maxvals'))';
        mdl.J_c = ((mdl.J_c')./(maxvals'))';
    end
    
    
    % compute discrete-time system matrices
    mdl.A = expm(mdl.A_c*mdl.dt);
    mdl.B = (mdl.A - eye(mdl.nx))*inv(mdl.A_c)*mdl.B_c;
    mdl.G = mdl.G_c;
    mdl.J = mdl.J_c;


end







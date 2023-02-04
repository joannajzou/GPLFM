function mdl = ss_MDOF(M,K,C,S_acc,S_p,dt)
% returns state-space matrices describing the model, measurements, and input

% INPUT:
%   M (ndof-by-ndof) = model mass matrix 
%   K (ndof-by-ndof) = model stiffness matrix 
%   C (ndof-by-ndof) = model damping matrix 
%   S_acc (na-by-ndof) = indicator matrix of measured accelerations
%   S_p (ndof-by-np) = indicator matrix of applied external loads
%   dt (float) = time step

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
    mdl.S_acc = S_acc;
    mdl.S_p = S_p;
    mdl.dt = dt;
    
    % modal frequencies
    eigenvals = eig(K,M);
    w = sqrt(eigenvals);
    f = w/(2*pi);
    
    % define dimensions
    mdl.ndof = size(M,1);                                            % number of dofs
    mdl.nm = 1;                                                      % number of forces
    mdl.nx = 2*mdl.ndof;                                             % number of states      
    mdl.ny = size(S_acc,1);                                          % number of measured states
    mdl.nye = mdl.ndof;                                              % number of outputs 
    mdl.np = size(S_p,2);                                            % number of external forces
    

    % compute continuous-time system matrices
    mdl.A_c = [zeros(mdl.ndof) eye(mdl.ndof);...
              -inv(mdl.M)*mdl.K -inv(mdl.M)*mdl.C];
      
    mdl.B_c = [zeros(mdl.ndof,mdl.np);...
              inv(mdl.M)*mdl.S_p];
    
    mdl.G_c = [-mdl.S_acc*inv(mdl.M)*mdl.K -mdl.S_acc*inv(mdl.M)*mdl.C];
    mdl.J_c = [mdl.S_acc*inv(mdl.M)*mdl.S_p];
    
    % for estimating all acceleration outputs
    mdl.G_e = [-inv(mdl.M)*mdl.K -inv(mdl.M)*mdl.C];    
    mdl.J_e = [inv(mdl.M)*mdl.S_p];
    
    
    % compute discrete-time system matrices
    mdl.A = expm(mdl.A_c*mdl.dt);
    mdl.B = (mdl.A - eye(mdl.nx))*inv(mdl.A_c)*mdl.B_c;
    mdl.G = mdl.G_c;
    mdl.J = mdl.J_c;


end







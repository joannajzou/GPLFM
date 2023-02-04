function [M,K,C] = model_ShearBldg(ndof, mass, k, dampratios, concmass)
% returns EOM matrices for a n-dof shear building model
% with mass proportional damping

% INPUT:
% ndof (int) = number of degrees of freedom (floors)
% mass (float) = mass at each floor
% k (float) = stiffness at each floor
% dampratios (nm-vector) [OPT] = damping ratios of the first nm modes
% concmass (float) [OPT] = concentrated mass added to top level 

% OUTPUT:
% M (ndof-by-ndof) = mass matrix
% K (ndof-by-ndof) = stiffness matrix
% C (ndof-by-ndof) = damping matrix


    % define EOM matrices
    M = mass*eye(ndof);
    K = 2*k*eye(ndof) + diag(-k*ones(1,ndof-1),1) + diag(-k*ones(1,ndof-1),-1);
    K(ndof,ndof) = k;
    
    % if defined, add concentrated mass
    if exist('concmass', 'var')
        M(end,end) = concmass;
    end
    
    % compute eigenvalue problem
    [Phi, Omega2] = eig(K,M);
    w = sqrt(diag(Omega2));         % natural frequencies (rad)
    f = w/(2*pi);                   % natural frequencies (Hz)
    zeta = 0.05*ones(size(f));      % damping ratio
    
    % if defined, change modal damping ratios
    if exist('dampratios', 'var')
        nm = length(dampratios);
        zeta(1:nm) = dampratios;
    end
    
    % compute damping matrix 
    Cg = diag(2*zeta.*w);            % damping matrix 
    C = M*Phi*Cg*Phi'*M;     

end





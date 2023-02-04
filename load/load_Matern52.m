function loads = load_Matern52(hp, dt, t_end, dofs)
% returns load struct for a harmonic force of the form f(t) = A*sin(2*pi*lamb*t)

% INPUT:
% hp (2-by-1) = Matern 5/2 hyperparameters [alpha, ls]
% dt (float) = time step (s)
% t_end (float) = end time of recording (s)
% dofs (list) = list of integer dofs where force is applied

% OUTPUT:
% load (struct) = contains load info

    loads = struct;

    nsteps = t_end/dt + 1;
    [knl, ~] = kernel_Matern52(hp, dt);

    % generate realization of GP with Matern kernel
    fp = zeros(1,nsteps); 
    z = zeros(knl.nz,1);
    for t = 2:nsteps 
    w = normrnd(0,knl.sigma_w);
    z_new = knl.F*z + knl.L*w;
    fp(t) = knl.H*z_new;
    z = z_new;
    end
    time = 0:dt:(nsteps-1)*dt;


    loads.type = "Matern52";
    loads.loc = dofs;
    loads.params = hp;
    loads.dt = dt;
    loads.time = time;
    loads.force = fp;


end
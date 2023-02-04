function loads = load_Harmonic(ampl, lamb, dt, t_end, dofs)
% returns load struct for a harmonic force of the form f(t) = A*sin(2*pi*f*t)

% INPUT:
% ampl (float) = amplitude of sinusoidal force
% f (float) = natural frequency
% dt (float) = time step (s)
% t_end (float) = end time of recording (s)
% dofs (list) = list of integer dofs where force is applied

% OUTPUT:
% load (struct) = contains load info

    loads = struct;

    nsteps = t_end/dt + 1;

    time = 0:dt:(nsteps-1)*dt;
    fp = zeros(1,nsteps);
    for t = 1:nsteps
        fp(t) = ampl*sin(2*pi()*lamb*(time(t)));     % sinusoidal force
    end


%     bound = 0.1;
%     range = 0:0.01:bound;
%     p = polyfit([0 bound/2 bound], [1 ampl/4+1 ampl+1], 2);
%     prof_pos = polyval(p, range);
%     prof_neg = flip(prof_pos);
%     
%     fp = zeros(1,nsteps);
%     for t = 1:length(range)
%         [time, fpt] = input_Harmonic(nsteps, dt, prof_pos(t), f+range(t)-bound);
%         [time, fpt2] = input_Harmonic(nsteps, dt, prof_neg(t), f+range(t));
%         fp = fp + fpt + fpt2;
%     end
    
    loads.type = "Harmonic";
    loads.loc = dofs;
    loads.params = [ampl, lamb];
    loads.dt = dt;
    loads.time = time;
    loads.force = fp;


end
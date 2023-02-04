function [t_n, x_n] = downsample(t, x, dt_n, t0_n)
% resamples a vector using new starting point and new time step,
% where dt_n >= dt

% INPUT:
% t (1-by-nt) = old time vector
% x (nx-by-nt) = old state vector (to be resampled)
% dt_n (scalar) = new time step
% t0_n (nx-by-1) [OPT] = new initial time

% OUTPUT:
% t_n (1-by-nt_n) = new time vector
% x_n (nx-by-nt_n) = new state vector

    
    nx = size(x,1);
    
    t_n = t(1):dt_n:t(end);
    nt_n = length(t_n);
    
    % interpolate
    x_n = zeros(nx, nt_n);
    for i = 1:nx
        x_n(i,:) = interp1(t, x(i,:), t_n);
    end
    
    % change initial time (truncates signal)
    if exist('t0_n','var')
        idx = find(t_n >= t0_n);
        t_n = t_n(idx);
        x_n = x_n(idx);
    end
    
end

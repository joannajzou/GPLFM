function z = hp_rescale(x,scale_old,scale_new)
% rescales vector according to new bounds

% INPUT:
% x (n-by-nx) = vector to be rescaled, where n=number of hps and nx=number of samples
% scale_old (n-by-2) [OPT] = lower and upper bounds of original scale of x
%              (if set equal to [], assumes x is normalized from 0 to 1)
% scale_new (n-by-2) [OPT] = lower and upper bounds of new scale of x
%              (if not provided, assumes x is to be normalized from 0 to 1)

% OUTPUT:
% z (n-by-1) = vector equivalent to x in the new scale


    n = size(x,1);
    
    % if new scale not provided, normalize
    if ~exist('scale_new','var')
        scale_new = [zeros(n,1), ones(n,1)];
    end

    % if old scale not provided, assume variable is normalized
    if isequal(scale_old,[])
        scale_old = [zeros(n,1), ones(n,1)];
    end

    % rescale
    z = zeros(size(x));
    for i = 1:n
        min = scale_old(i,1);
        max = scale_old(i,2);
        a = scale_new(i,1);
        b = scale_new(i,2);
        z(i,:) = (b-a)*(x(i,:) - min)/(max - min) + a;
    end

    
end
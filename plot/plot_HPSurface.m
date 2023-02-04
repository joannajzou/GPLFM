function [X,Y,Z] = plot_HPSurface(limits, n, objfunc)
% plots 3D mesh surface

% INPUT:
% limits (nhp-by-2) = array of lower and upper bounds of search domain
% n (int) = surface discretization
% objfunc (func) = objective function for HP tuning

% OUTPUT:
% X (n-by-n) = meshgrid values for dimension 1
% Y (n-by-n) = meshgrid values for dimension 2
% Z (n-by-n) = surface values from evaluating the objfunc


    range1 = linspace(limits(1,1), limits(1,2), n);
    range2 = linspace(limits(2,1), limits(2,2), n);

    [X,Y] = meshgrid(range1, range2);
    Z = zeros(n,n);
    
    for i = 1:n
        for j = 1:n
            disp([i,j]);
            val1 = X(i,j);
            val2 = Y(i,j);
            hpi = [val1; val2];
            hpi_norm = hp_rescale(hpi, [range1(1), range1(end); range2(1), range2(end)]);
            Z(i,j) = objfunc(hpi_norm);
        end
    end

    
end
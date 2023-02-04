function [x, fval, x_hist] = hp_fmincon(x0, objfun, A, b, Aeq, beq, ll, ul)
% wrapper function about fmincon 

% INPUT:
% x0 (nx-by-1) = initial state
% objfun (func) = objective function
% A, b, Aeq, beq, ll, ul = fmincon inputs

% OUTPUT:
% x (nx-by-1) = state at optimum
% fval (float) = optimum value 
% x_hist = gradient descent history


    x_hist = [];
    options = optimset('fmincon');
    options = optimset('OutputFcn', @outfun);
    options = optimset(options,'LargeScale','off');
    options = optimset(options,'Display','Iter');
    options = optimset(options,'MaxIter',200,'TolX',1e-12,'TolCon',1e-12,'TolFun',1e-12);
    
    [x, fval] = fmincon(objfun, x0, A, b, Aeq, beq, ll, ul, [], options);
    
    
    function stop = outfun(x,optimvalues,state)
        stop = false;
        if isequal(state,'iter')
            x_hist = [x_hist, x];
        end
    end

end
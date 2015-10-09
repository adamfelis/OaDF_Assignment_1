% Backtracking Linesearch
% Assumes that f and df are in the scope.

function steplength = BacktrackingLineSearch(xk, direction, steplength)
    %f (xk + ?pk) ? f (xk) + c?? f Tk pk
    p_low = 0.25;
    p_hi = 1 - p_low;
    c = 0.9;
    
    t1 = f(xk) + steplength * direction;
    t2 = f(xk) + c * steplength * df(xk)' * direction;
    
    while(f(xk + steplength * direction) > f(xk) + c * steplength * df(xk)' * direction)
        steplength = steplength * ((p_low + p_hi) / 2);
    end
   
end
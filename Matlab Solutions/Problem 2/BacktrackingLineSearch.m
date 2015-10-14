% Backtracking Linesearch

function steplength = BacktrackingLineSearch(xk, direction, steplength, f, df)
    p_low = 0.25;
    p_hi = 1 - p_low;
    c = 0.9;
    
    t0 = df(xk)';
    t01 = f(xk);
    t1 = f(xk) + steplength * direction;
    t2 = f(xk) + c * steplength * df(xk)' * direction;
    
    while(f(xk + steplength * direction) <= f(xk) + c * steplength * df(xk)' * direction)
        steplength = steplength * ((p_low + p_hi) / 2);
    end
   
end
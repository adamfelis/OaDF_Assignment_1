% Backtracking Linesearch

function steplength = BacktrackingLineSearch(xk, direction, steplength, f, df)
    p = 0.5;
    c = 0.2;
    
    while(f(xk + steplength * direction) > f(xk) + c * steplength * df(xk)' * direction)
        steplength = steplength * p;
    end
    
   
end
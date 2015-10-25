function [solution, information] = Newton(f, df, d2f, x_0, tolerance_for_Newton_algorithm)
    X_k                 = x_0; % Starting point
    max_iterations      = 100*length(X_k);
    %step_length = 1;

    % Initial iteration
    f_val = f(X_k);
    df_val = df(X_k);
    d2f_val = d2f(X_k);
    converged = (norm(df_val,'inf') <= tolerance_for_Newton_algorithm);
    iterations = 1;

    % Initial information data
    information.converged = converged;
    information.iterations = 1;
    information.X = X_k;
    information.F = f_val;
    information.dF = df_val;
    information.d2F = d2f_val;

    % Main loop
    while ~converged && (iterations < max_iterations)
        % ================================================
        
        %step_direction = -(inv(d2f_val))*df_val;
        step_direction = -d2f_val\df_val;
        step_length = 1;
        X_k_1 = X_k + step_length * step_direction;

        
        % The strong Wolfe Conditions line search
        c1 = 0.25;
        c2 = 1 - c1;
        rho = 0.9;
        while ( ...
            f(X_k_1) > f(X_k) + c1 * step_length * (df_val') * step_direction ...
            || ...
            abs((df(X_k_1)') * step_direction) < c2 * abs((df_val') * step_direction) ...
            )
            step_length = step_length * rho;
            X_k_1 = X_k + step_length * step_direction;
        end
        

        %step_length = BacktrackingLineSearch(X_k, step_direction, step_length, f, df);
        %X_k_1 = X_k + step_length * step_direction;

        f_val = f(X_k_1);
        df_val = df(X_k_1);
        d2f_val = d2f(X_k_1);
        converged = (norm(df_val,'inf') <= tolerance_for_Newton_algorithm);
        iterations = iterations+1;

        % Information data about the algorithm
        information.iterations = information.iterations+1;
        information.converged = converged;
        information.X = [information.X X_k_1];
        information.F = [information.F f_val];
        information.dF = [information.dF df_val];
        information.d2F = d2f_val;

        X_k = X_k_1;
    end
    solution = X_k_1;
end
    
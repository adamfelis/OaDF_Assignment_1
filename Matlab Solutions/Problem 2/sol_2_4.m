
% Newton algorithm with Newton direction
%   f = @(x1,x2) (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
%   df = {@(x1,x2) 2*x1 + 4*x1*(x1^2 + x2 - 11) + 2*x2^2 - 14}
%        {@(x1,x2) [2*x2 + 4*x2*(x2^2 + x1 - 7) + 2*x1^2 - 22}
%   They are defined in solution 2_0


%This initiates f, df, d2f
sol_2_0;

% Starting points
x_0 = [10 , 3 ;...
       -3 , 3 ;...
       0 , 0];
dimension = size(x_0);
amount_of_starting_points = dimension(1);

for i = 1 : 1 : amount_of_starting_points
    X_k                 = x_0(i, :); % Starting point
    max_iterations      = 100*length(X_k);
    converged           = false;
    stat.converged      = false; % converged
    stat.nfun           = 0; % number of function calls
    stat.iter           = 0; % number of iterations

    % Initial iteration
    f_val = f(X_k);
    df_val = df(X_k);
    d2f_val = d2f(X_k);
    converged = (norm(df_val,'inf') <= tolerance_for_Newton_algorithm);
    stat.nfun = 1;
    iterations = 1;

    % Store data for plotting
    stat.X = X_k;
    stat.F = f_val;
    stat.dF = df_val;
    stat.d2F = d2f_val;

    % Main loop of steepest descent
    while ~converged && (iterations < max_iterations)
        % Steepest descent step
        % ================================================
        step_direction = -(inv(d2f_val))*df_val;
        step_length = 1;
        X_k_1 = X_k + step_length * step_direction';


        % The strong Wolfe Conditions
        c1 = 0.25;
        c2 = 1 - c1;
        rho = 0.9;
        while ( ...
            f(X_k_1) > f(X_k) + c1 * step_length * (df_val') * step_direction ...
            || ...
            abs((df(X_k_1)') * step_direction) < c2 * abs((df_val') * step_direction) ...
            )
            step_length = step_length * rho;
            X_k_1 = X_k + step_length * step_direction';
        end


        % ================================================
        f_val = f(X_k_1);
        df_val = df(X_k_1);
        d2f_val = d2f(X_k_1);
        converged = (norm(df_val,'inf') <= tolerance_for_Newton_algorithm);

        iterations = iterations+1;
        stat.nfun = stat.nfun+1;

        % Store data for plotting
        stat.X = [stat.X X_k_1];
        stat.F = [stat.F f_val];
        stat.dF = [stat.dF df_val];
        stat.d2F = d2f_val;

        X_k = X_k_1;
    end

    % Prepare return data
    if ~converged
    X_k = [];
    end

    % If it converged
    stat.converged = converged;
    stat.iter = iterations;
    stat.ek = [];
    stat.gradient_norm = [];
    for i = 1:1:length(stat.X)
        stat.ek = [stat.ek norm(stat.X(i) - stat.X(length(stat.X)),2)];
        stat.gradient_norm = [stat.gradient_norm norm(stat.dF(i),'inf')];
    end
end

%plot(stat.X, stat.F);

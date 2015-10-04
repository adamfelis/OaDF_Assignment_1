function [x_out, converged, solution_info] = own_newton(x_in, fun, varargin)
addpath('immoptibox\');

% Initialize f, df, d2f, tolerance
sol_2_0;

%local parameters
converged = false;
operations_counter = 0;
operations_max_amount = 1500;
solution_vector_with_steps_values = [];
solution_vector_with_function_values = [];
%=================

x_old = x_in;
x_new = x_old;

if (norm(df,'inf') < tol)
    converged = true;
end

while(~converged && operations_counter < operations_max_amount)
    %newton_direction = -(d2f\df); 
    newton_direction = -df; 
    
    %[x_new, f_new, df_new, info] = linesearch(fun, x_old, f, df, newton_direction);%, 0, varargin{:});
    [x_new, info] = ucminf(fun, x_old);%, varargin{:});
    %options = optimoptions('fminunc','GradObj','on');
    %[x_new, f_new] = fminunc(fun, x_old, options);%, varargin{:});

    [f, df, d2f] = feval(fun, x_new);%, 0,0, varargin{:});
    solution_vector_with_steps_values = [solution_vector_with_steps_values x_new];
    solution_vector_with_function_values = [solution_vector_with_function_values f];
    if( (x_new - x_old)^2 < tol)
        disp('ERROR!!!!!!!!!!!!!!!!!');
        converged = false;
        break;
    end
    x_old = x_new;
    
    operations_counter = operations_counter + 1;
    converged = (norm(df,'inf') < tol);
end

x_out = x_new;


solution_info = table( transpose(solution_vector_with_steps_values), transpose(solution_vector_with_function_values));
end
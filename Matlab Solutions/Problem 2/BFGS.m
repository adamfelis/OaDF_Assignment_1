function [x_n, information] = BFGS( f, df, x_0, tolerance, max_amount_of_iterations, with_inversing_of_B )
% BFGS-based Quasi-Newton Algorithm
%   Input parameters:
%       f - function
%       df - derivative of function
%       x_0 - starting point
%       tolerance - if new solution is in the 'tolerance' neighbourhood of
%           previous solution algorithm will stop main loop execution (the best
%           solution is found)
%       max_amount_of_iterations - maximum amount of executed loops - if
%           reached, algorithm decides to stop execution without solution
%       with_inversing_of_B - if 'true' algorithm will inverse B matrix in
%           every loop; if 'false algorithm will compute H matrix
%   Output parameters:
%       x_n - final result
%       information - structure with some important facts about algorithm
%           execution
%           information.converged - boolean variable, which indicates whether 
%               algorithm converged or not
%           information.amount_of_iterations - amount of iterations used in
%               the algorithm
%           information.approximations - all aproximations of the result, 
%                collected through loop execution;

% Local variables:
iteration_counter   =        0;
B                   =        eye(length(x_0));
beta                =        1;
H                   =        beta * eye(length(x_0));
df_value            =        df(x_0);
converged           =        (norm(df_value, 'inf') < tolerance);
x_previous          =        x_0;
x_new               =        x_previous;
approximations      =        x_0;
%----------------------

% Main loop of the algorithm

while ~converged && (iteration_counter < max_amount_of_iterations)
    iteration_counter = iteration_counter + 1;
    
    % It works like: step = - (B)^(-1) * df_value
    if(with_inversing_of_B)
        step_direction = - (B \ df_value);
    else
        step_direction = - H * df_value;
    end
    step_length = 1;
    x_new = x_previous + step_length * step_direction;
    
    % The strong Wolfe Conditions
    eig(H)
    c1 = 10e-4;
    c2 = 0.9;%1 - c1;
    rho = 0.9;
    k = 0;
    while ( ...
        f(x_new) > f(x_previous) + c1 * step_length * (df_value') * step_direction ...
        || ...
        abs(step_direction' * df(x_new)) < c2 * abs(step_direction' * df_value) ...
        )
        step_length = step_length * rho;
        x_new = x_previous + step_length * step_direction;
        k = k+1;
    end
        disp(k);
    
    df_new_value = df(x_new);
    norm_value = norm(df_new_value, 'inf');
    converged = (norm_value <= tolerance);
    
    distance_in_arguments = x_new - x_previous;
    distance_in_derivative_values = df_new_value - df_value;
    
    if(with_inversing_of_B)
        B_times_distance = B * distance_in_arguments;
        B = B - (B_times_distance / dot(distance_in_arguments, B_times_distance)) * B_times_distance' + ...
            (distance_in_derivative_values / dot(distance_in_derivative_values, distance_in_arguments)) * distance_in_derivative_values';
    else
        
        if iteration_counter == 1
            beta = (distance_in_derivative_values' * distance_in_arguments) / (distance_in_derivative_values' * distance_in_derivative_values);
            H = beta * H;
        end
        
        ro = 1 / dot(distance_in_derivative_values, distance_in_arguments);
        E = ( eye(length(x_0)) - ro * distance_in_arguments * distance_in_derivative_values');
        H = E * H * E' + ro * (distance_in_arguments * distance_in_arguments');
    end
    x_previous = x_new;
    df_value = df_new_value;
    
    approximations = [ approximations x_new ];
end

information.converged = converged;
information.amount_of_iterations = iteration_counter;
information.approximations = approximations;
x_n = x_new;

if(information.converged)
    disp('Solution found for: ');
    disp(x_0);
    disp('f:');
    disp(f(x_n));
    disp('df:');
    disp(df(x_n));
else
    disp('Solution not found');
end

clear iteration_counter B beta H df_value converged x_previous x_new approximations;

end


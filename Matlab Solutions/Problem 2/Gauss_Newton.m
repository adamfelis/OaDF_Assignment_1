function [x_n, information] = Gauss_Newton(f, df, d2f, x_0, tolerance, max_amount_of_iterations)
% Gauss-Newton algorithm
%   Input arguments:
%       f - function
%       df - derivative of function calculated as J(x)' * r(x) 
%       d2f - second derivative of function calculated as J(x)' * J(x)
%       x_0 - starting point
%       tolerance - if new solution is in the 'tolerance' neighbourhood of
%           previous solution algorithm will stop main loop execution (the best
%           solution is found)
%       max_amount_of_iterations - maximum amount of executed loops - if
%           reached, algorithm decides to stop execution without solution
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


%% Local variables:
iteration_counter   =        0;
converged           =        (norm(df(x_0), 'inf') < tolerance);
x_previous          =        x_0;
x_new               =        x_previous;
approximations      =        x_0;

%% Main loop of the algorithm
while ~converged && (iteration_counter < max_amount_of_iterations)
    iteration_counter = iteration_counter + 1;
    
    
    %step_direction = - ((J'*J)^(-1)) * J' * r(x_previous);
    step_direction = - (d2f(x_previous)^(-1)) * df(x_previous);
    step_length = 1;
    x_new = x_previous + step_length * step_direction;
    
    % Armijo condition
    c1 = 0.25;
    rho = 0.9;
    while ( ...
        f(x_new) > f(x_previous) + c1 * step_length * (df(x_previous)') * step_direction ...
        )
        step_length = step_length * rho;
        x_new = x_previous + step_length * step_direction;
    end
    converged = (norm(df(x_new), 'inf') < tolerance);
    x_previous = x_new;
    approximations = [ approximations x_new ];
end

%% Final result
information.converged = converged;
information.amount_of_iterations = iteration_counter;
information.approximations = approximations;
x_n = x_new;

end


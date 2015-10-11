function [x_n, information] = Gauss_Newton(r, J, x_0, max_amount_of_iterations)
% Gauss-Newton algorithm
%   Input arguments:
%

%% Local variables:
iteration_counter   =        0;
converged           =        (norm(J' * r(x_0), 'inf') < tolerance);
x_previous          =        x_0;
x_new               =        x_previous;
approximations      =        x_0;

%% Main loop of the algorithm
while ~converged && (iteration_counter < max_amount_of_iterations)
    iteration_counter = iteration_counter + 1;
    
    
    step_direction = - ((J'*J)^(-1)) * J' * r(x_previous);
    step_length = 1;
    x_new = x_previous + step_length * step_direction;
    
    % Armijo condition
    c1 = 0.25;
    rho = 0.9;
    while ( ...
        f(x_new) > f(x_previous) + c1 * step_length * (r(x_previous)') * (J') * step_direction ...
        )
        step_length = step_length * rho;
        x_new = x_previous + step_length * step_direction;
    end
    converged = (norm(J' * r(x_0), 'inf') < tolerance);
    x_previous = x_new;
end


end


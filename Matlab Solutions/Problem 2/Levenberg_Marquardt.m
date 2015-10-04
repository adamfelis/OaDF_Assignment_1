function [x_n, information] = Levenberg_Marquardt( f, df, d2f, x_0, r, tolerance, max_amount_of_iterations )
% Levenberg_Marquardt Algorithm
%   Input parameters:
%       f - function
%       df - derivative of function
%       d2f - second derivative of function
%       x_0 - starting point
%       r - value used to set up mi_0 parameter, when >=1 then the initial
%           matrix is guaranteed to be positive defined; default: 10e-6
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

% Local variables:
iteration_counter               =   0;
x                               =   x_0;
mi_0                            =   r * norm(d2f(x), 'inf');
mi                              =   mi_0;
ni                              =   2;
size_of_hessian                 =   length(d2f(x));
A                               =   d2f(x) + mi * eye(size_of_hessian);
found                           =   false;
m_k                             =   @(x_k, h, B) f(x_k) + (df(x_k)') * (h) + 0.5 * h' * B * h;
delta                           =   0;
approximations                  =   x_0;
%----------------------

while ( ~found && iteration_counter <=  max_amount_of_iterations )
    
    A = d2f(x) + mi * eye(size_of_hessian);
    is_matrix_A_positive_defined = is_matrix_positive_defined(A); 
    while( ~is_matrix_A_positive_defined )
        mi = 2 * mi;
        A = d2f(x) + mi * eye(size_of_hessian);
        is_matrix_A_positive_defined = is_matrix_positive_defined(A); 
    end
    
    % Solving equation A * h_dn = - df ; A = R'*R - cholesky factorization
    R = chol(A);
    h_dn = R\(R'\(-df(x)));
    ro = (f(x) - f(x + h_dn)) / (m_k(x, [0;0], d2f(x)) - m_k(x, h_dn, d2f(x)));
    if ( ro > delta )
        x = x + h_dn;
        approximations = [ approximations x ];
        ni = 2;
        mi = mi * max(1/3, 1 - (2*ro - 1)^3);
    else
        mi = mi * ni;
        ni = 2 * ni;
    end
    iteration_counter = iteration_counter + 1;
    
    inf_norm_value = norm(df(x), 'inf');
    found = (inf_norm_value <= tolerance);
end

information.converged = found;
information.amount_of_iterations = iteration_counter;
information.approximations = approximations;

x_n = x;

end
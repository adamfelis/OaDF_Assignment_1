% 2.7: Levenberg-Marquardt algorithm
%   Implement the Levenberg-Marquardt algorithm. You should motivate your algorithmic
%   choices and explain the algorithm. Test this algorithm for various starting points.
%   Use the contour plot to demonstrate the convergence sequence for this algorithm.
%   Make tables and plots of the convergence rate. Is the convergence rate as expected.

sol_2_0;

% Script parameters:
% Rows of x_0 are starting points
x_0 = [10 , 3 ;...
       -3 , 3 ;...
        0 , 0];
dimension = size(x_0);
amount_of_starting_points = dimension(1);
%--------------------

% Global variables:
final_results_of_2_7 = struct([]);
%-------------------

[x_n, information] = Levenberg_Marquardt(f, df, d2f, x_0(2, :)', 0.01, tolerance_for_Levenberg_Marquardt_algorithm, max_amount_of_iterations);

%clearvars -except final_results_of_2_7;
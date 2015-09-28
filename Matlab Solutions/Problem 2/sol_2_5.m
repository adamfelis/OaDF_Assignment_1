% 2.5: Quasi-Newton algorithm
%   Implement a BFGS-based Quasi-Newton algorithm with a line search algorithm
%   of your choice. You should motivate your algorithmic choices and explain the algorithm.
%   Test this algorithm for various starting points. Use the contour plot to
%   demonstrate the convergence sequence for this algorithm. Make tables and plots of
%   the convergence rate. Is the convergence rate as expected.

sol_2_1;

% Script parameters:
% Rows of x_0 are starting points
x_0 = [10 , 3 ;...
       -3 , 3 ;...
        0 , 0];
dimension = size(x_0);
amount_of_starting_points = dimension(1);
%--------------------

% Global variables:
final_results_of_2_5 = struct([]);
%-------------------

for i = 1 : 1 : amount_of_starting_points;
    [x_n, information] = BFGS(f, df, x_0(i, :)', tolerance_for_BFGS_algorithm, max_amount_of_iterations, true);
    final_results_of_2_5(i).information = information;
    final_results_of_2_5(i).solution = x_n;
    if( ~information.converged)
        continue;
    end
    
    figure(i+1)
    subplot(1,2,1);
    [c,h] = contour(X,Y,Composition_matrix_for_z_axis,v,'linewidth',2);
    colorbar;
    axis image;
    xlabel('x_1','Fontsize',14);
    ylabel('x_2','Fontsize',14);
    hold on;
    plot(information.approximations(1,:), information.approximations(2,:));
    
    subplot(1,2,2);
    final_results_of_2_5(i).convergence_rates = (plot_convergence(information.approximations, x_n))';
    
end

clearvars -except final_results_of_2_5;
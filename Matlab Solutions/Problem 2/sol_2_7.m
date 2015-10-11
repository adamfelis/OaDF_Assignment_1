% 2.7: Levenberg-Marquardt algorithm
%   Implement the Levenberg-Marquardt algorithm. You should motivate your algorithmic
%   choices and explain the algorithm. Test this algorithm for various starting points.
%   Use the contour plot to demonstrate the convergence sequence for this algorithm.
%   Make tables and plots of the convergence rate. Is the convergence rate as expected.

sol_2_1;

%% Script parameters:
% Rows of x_0 are starting points
x_0 = [10 , 3 ;...
       -3 , 3 ;...
        0 , 0];
dimension = size(x_0);
amount_of_starting_points = dimension(1);
%--------------------

%% Global variables:
final_results_of_2_7 = struct([]);
%-------------------

%% Main Loop
for i = 1 : 1 : amount_of_starting_points
    [x_n, information] = Levenberg_Marquardt(f, df, d2f, x_0(i, :)', 10e-6, tolerance_for_Levenberg_Marquardt_algorithm, max_amount_of_iterations);
    final_results_of_2_7(i).information = information;
    final_results_of_2_7(i).solution = x_n;
    if( ~information.converged)
        continue;
    end
    
    figure(i+1)
    subplot(2,2,[1 3]);
    [c,h] = contour(X,Y,Composition_matrix_for_z_axis,v,'linewidth',2);
    colorbar;
    axis image;
    xlabel('x_1','Fontsize',14);
    ylabel('x_2','Fontsize',14);
    hold on;
    plot(information.approximations(1,:), information.approximations(2,:));
    plot(information.approximations(1,:), information.approximations(2,:), 'xr');
    
    subplot(2,2,2);
    final_results_of_2_7(i).convergence_rates = (plot_convergence(information.approximations, x_n))';
    subplot(2,2,4);
    plot_convergence_rate(information.approximations, x_n);
end

clearvars -except final_results_of_2_7;
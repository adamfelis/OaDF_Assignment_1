
% Steepest descent algorithm with steepest descent direction Bk-df with Bk = I
%   f = @(x1,x2) (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
%   df = {@(x1,x2) 2*x1 + 4*x1*(x1^2 + x2 - 11) + 2*x2^2 - 14}
%        {@(x1,x2) [2*x2 + 4*x2*(x2^2 + x1 - 7) + 2*x1^2 - 22}
%   They are defined in solution 2_0


%This initiates f, df, d2f
%sol_2_0;
sol_2_1;

% Starting points
x_0 = [10 , 3 ;...
       -3 , 3 ;...
       0 , 0];
dimension = size(x_0);
amount_of_starting_points = dimension(1);


for i = 1 : 1 : amount_of_starting_points
    [solution, information] = SteepestDescent(f, df, x_0(i,:)', tolerance_for_SDD_algorithm);

    % It didnt converge
    if ~information.converged
        continue;
    end

    % If it converged
    %{
    stat.converged = converged;
    stat.iter = iterations;
    stat.ek = [];
    stat.gradient_norm = [];
    for i = 1:1:length(stat.X)
        stat.ek = [stat.ek norm(stat.X(i) - stat.X(length(stat.X)),2)];
        stat.gradient_norm = [stat.gradient_norm norm(stat.dF(i),'inf')];
    end
    %}
    
    
    % PLOT
    figure(i+1)
    subplot(1,2,1);
    [c,h] = contour(X,Y,Composition_matrix_for_z_axis,v,'linewidth',2);
    colorbar;
    axis image;
    title(strcat('Starting point = ', mat2str(x_0(i,:))));
    xlabel('x_1','Fontsize',14);
    ylabel('x_2','Fontsize',14);
    hold on;
    test = information.X;
    plot(information.X(1,:), information.X(2,:));

    subplot(1,2,2);
    convergence_rates = (plot_convergence(information.X, solution))';
    % ===============================================================
    
end

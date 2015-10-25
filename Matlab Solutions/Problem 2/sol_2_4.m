
% Newton algorithm with Newton direction
%   f = @(x1,x2) (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
%   df = {@(x1,x2) 2*x1 + 4*x1*(x1^2 + x2 - 11) + 2*x2^2 - 14}
%        {@(x1,x2) [2*x2 + 4*x2*(x2^2 + x1 - 7) + 2*x1^2 - 22}
%   They are defined in solution 2_0


% NOTES: Currently Newton doesnt converge with starting point 0,0 with
% Wolfe linesearch
% When using Backtracking linesearch, one of the graphs are are wierd,
% check it.


%This initiates f, df, d2f
sol_2_1;
%sol_2_0;

% Starting points
x_0 = [10 , 3 ;...
       -3 , 3 ;...
       0 , 0];
dimension = size(x_0);
amount_of_starting_points = dimension(1);


for i = 1 : 1 : amount_of_starting_points
    [solution, information] = Newton(f, df, d2f, x_0(i,:)', tolerance_for_Newton_algorithm);

    % It didnt converge
    if ~information.converged
        continue;
    end

    % If it converged
    %{
    stat.ek = [];
    stat.gradient_norm = [];
    for ii = 1:1:length(information.X)
        stat.ek = [stat.ek norm(information.X(ii) - information.X(length(information.X)),2)];
        stat.gradient_norm = [stat.gradient_norm norm(information.dF(ii),'inf')];
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
    %convergence_rates = (plot_convergence(information.X, solution))';
    e = plot_convergence_rate(information.X, solution);
    solution_table = prepare_table(e, information.X, f, df, solution);
    % ===============================================================
end

clearvars
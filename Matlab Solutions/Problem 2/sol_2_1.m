% 2.1: Contour plot and stationary points
%   Make a contour plot of f(x) and locate all local minimizers, all local
%   maximizers, and all saddle points.

sol_2_0;

% Script parameters:
boundary_for_first_parameter.start      =    -5;
boundary_for_first_parameter.end        =    max(boundary_for_first_parameter.start, 5);
boundary_for_second_parameter.start     =    -5;
boundary_for_second_parameter.end       =    max(boundary_for_second_parameter.start, 5);

density_for_parameter                   =    1000;

step_for_first_parameter                =    (boundary_for_first_parameter.end - boundary_for_first_parameter.start) / density_for_parameter;
step_for_second_parameter               =    (boundary_for_second_parameter.end - boundary_for_second_parameter.start) / density_for_parameter;
% ------------------

% Global variables:
final_results_of_2_1 = struct();
%-------------------

% Solution:
x1 = boundary_for_first_parameter.start + step_for_first_parameter : step_for_first_parameter : boundary_for_first_parameter.end;
x2 = boundary_for_second_parameter.start + step_for_second_parameter : step_for_second_parameter : boundary_for_second_parameter.end;

[X,Y] = meshgrid(x1,x2);

Composition_matrix_for_z_axis = zeros(density_for_parameter, density_for_parameter);
for i = 1 : density_for_parameter
    for j = 1 : density_for_parameter
        Composition_matrix_for_z_axis(i,j) = f([X(i,j), Y(i,j)]);
    end
end

% Uncomment this part if You would like to see the figure with the function


figure(1)
v = [0:2:10 10:10:100 100:20:200];
[c,h] = contour(X,Y,Composition_matrix_for_z_axis,v,'linewidth',2);
colorbar;
axis image;
xlabel('x_1','Fontsize',14);
ylabel('x_2','Fontsize',14);

% Findong all local minimizers, maximizers and saddle points:


% Minimizers:
sol_1 = fsolve( df, [10 10], options);
sol_2 = fsolve( df, [-10 10], options);
sol_3 = fsolve( df, [-10 -10], options);
sol_4 = fsolve( df, [10 -10], options);
% Maximizers:
sol_5 = fsolve( df, [0 0], options);
% Saddle points:
sol_6 = fsolve( df, [-3 0], options);
sol_7 = fsolve( df, [3 0], options);
sol_8 = fsolve( df, [0 3], options);

final_results_of_2_1.minimizers = [sol_1 ; sol_2; sol_3; sol_4];
final_results_of_2_1.maximizers = sol_5;
final_results_of_2_1.saddle_points = [sol_6 ; sol_7; sol_8];

clearvars -except final_results_of_2_1 f df d2f options tolerance_for_BFGS_algorithm max_amount_of_iterations X Y v Composition_matrix_for_z_axis;
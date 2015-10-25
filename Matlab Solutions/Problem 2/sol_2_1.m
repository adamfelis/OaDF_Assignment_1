% 2.1: Contour plot and stationary points
%   Make a contour plot of f(x) and locate all local minimizers, all local
%   maximizers, and all saddle points.

sol_2_0;

%% Script parameters:
boundary_for_first_parameter.start      =    -5;
boundary_for_first_parameter.end        =    max(boundary_for_first_parameter.start, 5);
boundary_for_second_parameter.start     =    -5;
boundary_for_second_parameter.end       =    max(boundary_for_second_parameter.start, 5);

density_for_parameter                   =    100;

step_for_first_parameter                =    (boundary_for_first_parameter.end - boundary_for_first_parameter.start) / density_for_parameter;
step_for_second_parameter               =    (boundary_for_second_parameter.end - boundary_for_second_parameter.start) / density_for_parameter;
% ------------------

%% Global variables:
final_results_of_2_1 = struct();
%-------------------

%% Function plot:
x1 = boundary_for_first_parameter.start + step_for_first_parameter : step_for_first_parameter : boundary_for_first_parameter.end;
x2 = boundary_for_second_parameter.start + step_for_second_parameter : step_for_second_parameter : boundary_for_second_parameter.end;

[X,Y] = meshgrid(x1,x2);

Composition_matrix_for_z_axis = zeros(density_for_parameter, density_for_parameter);
for i = 1 : density_for_parameter
    for j = 1 : density_for_parameter
        Composition_matrix_for_z_axis(i,j) = f([X(i,j), Y(i,j)]);
    end
end


figure(1)
v = [0:2:10 10:10:100 100:20:200];
[c,h] = contour(X,Y,Composition_matrix_for_z_axis,v,'linewidth',2);
colorbar;
axis image;
xlabel('x_1','Fontsize',14);
ylabel('x_2','Fontsize',14);

%% Finding all local minimizers, maximizers and saddle points:

syms a b;
f_ab = (a^2 + b - 11)^2 + (a + b^2 - 7)^2;
df_da = diff(f_ab, a);
df_db = diff(f_ab, b);
[y1, y2] = solve(df_da == 0, df_db == 0, a, b);
minimizers = real([double(y1) double(y2)]);
hold on;
plot(minimizers(:,1), minimizers(:,2), '.r', 'MarkerSize', 25);

%% Result
final_results_of_2_1.minimizers = minimizers;

clearvars -except final_results_of_2_1 f df d2f options tolerance_for_BFGS_algorithm tolerance_for_Levenberg_Marquardt_algorithm max_amount_of_iterations X Y v Composition_matrix_for_z_axis tolerance_for_Newton_algorithm tolerance_for_SDD_algorithm tolerance_for_Gauss_Newton_algorithm;

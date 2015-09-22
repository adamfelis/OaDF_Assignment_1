% 2.1: Contour plot and stationary points
%   Make a contour plot of f(x) and locate all local minimizers, all local
%   maximizers, and all saddle points.

sol_2_0;

% Script parameters:
boundary_for_first_parameter.start = -5;
boundary_for_first_parameter.end = max(boundary_for_first_parameter.start, 5);
boundary_for_second_parameter.start = -5;
boundary_for_second_parameter.end = max(boundary_for_second_parameter.start, 5);

density_for_parameter = 2000;

step_for_first_parameter = (boundary_for_first_parameter.end - boundary_for_first_parameter.start) / density_for_parameter;
step_for_second_parameter = (boundary_for_second_parameter.end - boundary_for_second_parameter.start) / density_for_parameter;
% ------------------

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

figure(1)
v = [0:2:10 10:10:100 100:20:200];
[c,h] = contour(X,Y,Composition_matrix_for_z_axis,v,'linewidth',2);
colorbar;
axis image;
xlabel('x_1','Fontsize',14);
ylabel('x_2','Fontsize',14);
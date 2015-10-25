% 2.2: Gradient and Hessian
%   Derive the gradient and the Hessian for f(x). Verify that rf(x) = 0 for all stationary
%   points in 2.1. What property does the Hessian have for a) the local minimizers, b)
%   the local maximizers, and c) the saddle points. Verify this numerically.

sol_2_1;

%% Script parameters:
dfs = zeros(length(final_results_of_2_1.minimizers),2);
eigenvalues = zeros(length(final_results_of_2_1.minimizers),2);

%% Main loop:

for i = 1 : 1 : length(final_results_of_2_1.minimizers)
    sol = df(final_results_of_2_1.minimizers(i,:));
    dfs(i,1) = sol(1);
    dfs(i,2) = sol(2);
    
    values = eig(d2f(final_results_of_2_1.minimizers(i,:)));
    eigenvalues(i,1) = values(1);
    eigenvalues(i,2) = values(2);
end

%% Result
final_results_of_2_2.dfs = dfs;
final_results_of_2_2.eigenvalues = eigenvalues;

clearvars -except final_results_of_2_2
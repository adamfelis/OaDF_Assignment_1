% 2.2: Gradient and Hessian
%   Derive the gradient and the Hessian for f(x). Verify that rf(x) = 0 for all stationary
%   points in 2.1. What property does the Hessian have for a) the local minimizers, b)
%   the local maximizers, and c) the saddle points. Verify this numerically.

sol_2_1;

df_1 = df(final_results_of_2_1.minimizers(1).x);
df_2 = df(final_results_of_2_1.minimizers(2).x);
df_3 = df(final_results_of_2_1.minimizers(3).x);
df_4 = df(final_results_of_2_1.minimizers(4).x);

df_5 = df(final_results_of_2_1.maximizers(1).x);

df_6 = df(final_results_of_2_1.saddle_points(1).x);
df_7 = df(final_results_of_2_1.saddle_points(2).x);
df_8 = df(final_results_of_2_1.saddle_points(3).x);

Eigenvalues_Of_The_Hessian_For_The_First_Minimizer = eig(d2f(final_results_of_2_1.minimizers(1).x));
Eigenvalues_Of_The_Hessian_For_The_Second_Minimizer = eig(d2f(final_results_of_2_1.minimizers(2).x));
Eigenvalues_Of_The_Hessian_For_The_Third_Minimizer = eig(d2f(final_results_of_2_1.minimizers(3).x));
Eigenvalues_Of_The_Hessian_For_The_Fourth_Minimizer = eig(d2f(final_results_of_2_1.minimizers(4).x));


Eigenvalues_Of_The_Hessian_For_The_Maximizer = eig(d2f(final_results_of_2_1.maximizers(1).x));


Eigenvalues_Of_The_Hessian_For_The_First_Saddle_Point = eig(d2f(final_results_of_2_1.saddle_points(1).x));
Eigenvalues_Of_The_Hessian_For_The_Second_Saddle_Point = eig(d2f(final_results_of_2_1.saddle_points(2).x));
Eigenvalues_Of_The_Hessian_For_The_Third_Saddle_Point = eig(d2f(final_results_of_2_1.saddle_points(3).x));
% 2.2: Gradient and Hessian
%   Derive the gradient and the Hessian for f(x). Verify that rf(x) = 0 for all stationary
%   points in 2.1. What property does the Hessian have for a) the local minimizers, b)
%   the local maximizers, and c) the saddle points. Verify this numerically.

sol_2_1;

Eigenvalues_Of_The_Hessian_For_The_First_Minimizer = eig(d2f(sol_1));
Eigenvalues_Of_The_Hessian_For_The_Second_Minimizer = eig(d2f(sol_2));
Eigenvalues_Of_The_Hessian_For_The_Third_Minimizer = eig(d2f(sol_3));
Eigenvalues_Of_The_Hessian_For_The_Fourth_Minimizer = eig(d2f(sol_4));


Eigenvalues_Of_The_Hessian_For_The_Maximizer = eig(d2f(sol_5));


Eigenvalues_Of_The_Hessian_For_The_First_Saddle_Point = eig(d2f(sol_6));
Eigenvalues_Of_The_Hessian_For_The_Second_Saddle_Point = eig(d2f(sol_7));
Eigenvalues_Of_The_Hessian_For_The_Third_Saddle_Point = eig(d2f(sol_8));
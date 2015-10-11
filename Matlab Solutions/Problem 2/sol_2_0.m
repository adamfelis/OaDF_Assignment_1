% Default starting script for Problem 2

clear all;
close all;

f   = @(X)  ( (X(1).^2 + X(2) - 11).^2 + (X(1) + X(2).^2 - 7).^2 ); 
df  = @(X)  [ (4.*X(1).^3 + 4.*X(1).*X(2) - 42.*X(1) + 2.*X(2).^2 - 14); ...
              (2.*X(1).^2 - 26.*X(2) - 22 + 4.*X(1).*X(2) + 4.*X(2).^3) ];
d2f = @(X)  [ (12.*X(1).^2 + 4.*X(2) -42)   ,   (4.*X(1) + 4.*X(2)) ; ...
              (4.*X(1) + 4.*X(2))           ,   (12.*X(2)^2 + 4.*X(1) - 26)];
          
options = optimset('Display','off');


tolerance_for_SDD_algorithm = 1.0e-10;
tolerance_for_Newton_algorithm = 1.0e-10;
tolerance_for_BFGS_algorithm = 1e-10;
tolerance_for_Levenberg_Marquardt_algorithm = 1e-10;
tolerance_for_Gauss_Newton_algorithm = 1e-10;

% Modify this variable to modify maximum amount 
max_amount_of_iterations = 1000;
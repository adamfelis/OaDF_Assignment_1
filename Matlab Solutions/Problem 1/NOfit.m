function [ x_star, r_star ] = NOfit( t,y,n )
% NOfit function computes least squares fit with the model:
%
%           M(x,scalar_t) = x_1 + x_2*sin(w*scalar_t) + x_3*cos(w*scalar_t) 
%           + x_4*sin(2w*scalar_t) + x_5*cos(2w*scalar_t) + ...
%
% Input parameters:
%       - t - vector of data, which are known to the model; contains scalar
%           values which indicates time measurement. Further details are 
%           described below,
%       - y - vector of data, measured for each time moment from the vector
%           t (it means that measurement y(i) corresponds to t(i) moment in
%           time),
%       - n - value, which describes precision of the model (only first n
%           values from the sum in the model M are taken into account).
%
% The model M(x,scalar_t) is computed for each element of vector t 
% specified as input parameter. This vector indicates points in time, when 
% measurement y was obtained. For instance t = [1 4 5 7], means that y(1) 
% is a measurement after 1 hour, y(2) is measurement after 4 hours, etc. It
% means that the length of vector t is equal to the length of vector y. In 
% general elements of t not always are scalar and describe time. Those 
% could be also vectors which contain additional information (ex. element 
% of t could be a vector of two elements - one represents time and second
% one weight of the patient or his height in the case of health
% measurements on a big group of people).
%
% Assumptions:
%       - elements of t represents time (are scalars) and t(i-1) < t(i)
%       - length(y) = length(t)
%       - 3 <= n
%       - mod(n,2) = 1
%
% Function NOfit calculates so called residuals r = r(x) = Ax - y, where A 
% is nxn matrix (m = length(t)): 
%                   |1  sin(w*t(1))  cos(w*t(1))  sin(2w*t(1)) ...  |
%                   |1  sin(w*t(2))  cos(w*t(2))  sin(2w*t(2)) ...  |
%                   |1       .             .            .      ...  |
%               A = |1       .             .            .      ...  |
%                   |1       .             .            .      ...  |
%                   |1  sin(w*t(m))  cos(w*t(m))  sin(2w*t(m)) ...  |
% and minimizes function f(x) = 0.5 * norm(r(x),2)^2 = 0.5*norm(Ax - y,2)^2.
% As a result of minimization vector x_star is obtained.
% 
% Output parameters:
%       - x_star - n-length vector which minimizes f(x) function described
%           above (this vector provides the best model fitting),
%       - r_star - n-length vector of residuals computed for x_star
%           minimization vector.

%% Check whether input parameters fulfill assumptions:

if length(t) ~= length(y)
    disp('error - t and y have different lengths');
    x_star = [];
    r_star = [];
    return;
end
if n < 3
    disp('error - wrong n parameter specified - n < 3');
    x_star = [];
    r_star = [];
    return;
end
if mod(n,2) == 0
    disp('error - wrong n parameter specified - 2|n');
    x_star = [];
    r_star = [];
    return;
end

for i = 2: 1: length(t)
   if t(i) <= t(i - 1)  
       disp('error - t is not increasing');
       x_star = [];
       r_star = [];
       return;
   end
end

%% Local variables:
w                   =    (2 * pi) / 24;
m                   =    length(t);
A                   =    ones(m,n);

%% Main program:

% That loops fill matrix A by values specified in function description. 
% First column is not taken under account due to the fact that in the 
% beginning A is a matrix created with ones() and the first column is 
% already fulfilled correctly.
for i = 1 : 1 : m
    for j = 2 : 1 : n
        if mod(j, 2) == 0
            A(i,j) = sin((j/2)*w*t(i));
        else
            A(i,j) = cos(((j-1)/2)*w*t(i));
        end
    end
end

% Residual function is specified as a difference between a model (Ax) and
% measured values (y).
r = @(x) (A * x - y);

% Objective function  f could be written as follows: f(x) =
% 0.5*norm(r(x),2)^2. To calculate the minimizer Jacobian df/dx is 
% calculated:
%               df/dx = A'*(Ax-y).
% To obtain a minimizer we have to calulate solution of the equation:
%               df/dx = 0.
% Also: 
%           A'Ax -A'y = 0 -> A'Ax = A'y -> x = (A'A) \ (A')y

x_star = (A'*A) \ (A')*y;

% r_star is a residual vector for a minimum of the objective function.
r_star = r(x_star);

end
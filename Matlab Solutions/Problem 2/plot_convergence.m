function convergence_rates = plot_convergence( approximations, result )
% Method plots convergence diagram of all approximated values gathered, 
% through the execution of the approximation algorithm.
%   Input parameters:
%       approximations - approximated values
%       result - final result
%   Output parameters : 
%       convergence_rates - calculated rates of convergence for all
%           approximated values

% Local variables:
approximations_size =   size(approximations);
iterations_amount   =   approximations_size(2);
iterations          =   1:1:iterations_amount;
convergence_rates   =   zeros(1, iterations_amount); 
%-----------------

for i = 1 : 1 : iterations_amount
    convergence_rates(i) = norm(approximations(:,i) - result, 'inf');
end

plot(iterations, convergence_rates);

clear iterations_amount iterations;

end


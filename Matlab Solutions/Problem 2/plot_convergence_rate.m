function plot_convergence_rate(approximations, result)

% Local variables:
approximations_size =   size(approximations);
iterations_amount   =   approximations_size(2);
result_vector = ones(2, iterations_amount);
result_vector(1, : ) = result_vector(1, : ) * result(1);
result_vector(2, : ) = result_vector(2, : ) * result(2);
e = abs(approximations - result_vector);
log_e = log(e);
%-----------------
X = log_e(:,1:length(log_e)-1);
Y = log_e(:,2:length(log_e));
axis image;
plot(X(1,:), Y(2,:)); 
hold on;
plot(X(1,:), Y(2,:), 'xr'); 

end


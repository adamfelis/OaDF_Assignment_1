function table = prepare_table( e, approximations, f, df, sol)

amount_of_iterations = length(e);
amount_of_last_iterations = min(6, amount_of_iterations);

table = zeros(amount_of_last_iterations, 3);

for i = 1 : 1 : amount_of_last_iterations
    table(i,1) = norm(df(approximations(:,amount_of_iterations - (amount_of_last_iterations - i))),'inf');
   table(i,2) = norm(approximations(:,amount_of_iterations - (amount_of_last_iterations - i)) - sol,2);
   table(i,3) = abs(f(approximations(:,amount_of_iterations - (amount_of_last_iterations - i))) - f(sol));
end


end


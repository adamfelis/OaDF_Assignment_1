function is_positive = is_matrix_positive_defined( A )
% is_matrix_positive_defined checks whether matrix is positive defined
%   Input parameters:
%       A - matrix to check
%   Output parameters:
%       is_positive - boolean answer whether matrix A is positive defined

% Local variables:
eigenvalues             =   eig(A);
%-------------------

is_positive = true;

for i = 1 : 1 : length(eigenvalues)
    if(eigenvalues(i) <= 0)
        is_positive = false;
        break;
    end
end


end


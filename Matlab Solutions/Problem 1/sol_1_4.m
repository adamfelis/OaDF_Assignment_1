sol_1_0;

n = 13;

[x_star, r_star] = NOfit(t,Y,n);

m = 24;
w = (2 * pi) / 24;
A = ones(m,n);

for i = 1 : 1 : m
    for j = 2 : 1 : n
        if mod(j, 2) == 0
            A(i,j) = sin((j/2)*w*t(i));
        else
            A(i,j) = cos(((j-1)/2)*w*t(i));
        end
    end
end

cov_x = ((norm(r_star, 2) / sqrt(m-n))^2) * (A'*A)^(-1);

%% Initialization
sol_1_0;
n = 3;
t_continuous = linspace(t(1), t(length(t)), 120);
w = (2*pi)/24;

%% Main program
[x_star, r_star] = NOfit(t,Y,n);

eps = 0.001;
if norm(r_star, 2) < 292.558 - eps || norm(r_star, 2) > 292.558 + eps
    disp('Wrong result obtained by NOfit - second norm of r_star should be 292.558');
    return;
end

disp('Found solution: ');
disp(x_star);

M = @(x,t) (x(1) + x(2) * sin(w*t) + x(3) * cos(w*t));

figure(1);
set(gcf,'numbertitle','off','name','Least squares fit M(x*,t)') 
plot(t_continuous, M(x_star,t_continuous));
hold on;
plot(t,Y,'or');

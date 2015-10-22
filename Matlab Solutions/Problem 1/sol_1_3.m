sol_1_0;
n = 7;
[x_star, r_star] = NOfit(t,Y,n);

m = length(r_star);
n_plus = 0;
n_minus = 0;
runs = 0;
current_run = 'negative';


% This calculates the values needed to test for random signs
for i = 1 : 1 : m
    if strcmp(current_run, 'negative')
        if r_star(i) > 0
            n_plus = n_plus + 1;
            runs = runs + 1;
            current_run = 'positive';
            continue;
        end
        if r_star(i) < 0
            n_minus = n_minus + 1;
            continue;
        end
    end
    
    if strcmp(current_run, 'positive')
        if r_star(i) > 0
            n_plus = n_plus + 1;
            continue;
        end
        if r_star(i) < 0
            n_minus = n_minus + 1;
            runs = runs + 1;
            current_run = 'negative';
            continue;
        end
    end
end

mean = (2*n_plus*n_minus / m) + 1;  
deviation = ((mean-1)*(mean-2)) / (m-1);
test_for_random_signs_value = (abs(runs - mean)) / (sqrt(deviation));
random_likely = false;

if test_for_random_signs_value <= 1.96
    random_likely = true;
end

% This calculates the data needed to test for correlation
autocorrelation = 0;
first_part_of_trend_threshhold = 1 / (sqrt(m-1));
sum_part_of_threshhold = 0;
trend_treshhold = 0;
trend_likely = false; % There should be no trends in the residuals.

for i = 1 : 1 : m-1
    autocorrelation = autocorrelation + r_star(i)*r_star(i+1);
end

for i = 1 : 1 : m
    sum_part_of_threshhold = sum_part_of_threshhold + (r_star(i)^2);
end

trend_treshhold = first_part_of_trend_threshhold * sum_part_of_threshhold;

if abs(autocorrelation) > trend_treshhold
    trend_likely = true;
end

% Assume zero means, check this out
% Kig på de n værdier som opfylder at der ikke er trends og at der ikke er
% høj tilfældighed. Så kig på residual værdierne i dem og lav graf over dem
% så man kan se hvilken n der har residuals der reducerer mest.

y = 1;
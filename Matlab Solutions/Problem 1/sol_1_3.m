sol_1_0;
[dimension1,dimension2] = size(all_possible_n);
figure(1); % Residuals as function of n
figure(2); % Random signs as function of n
figure(3); % Trends as function of n
figure(4); % Residuals where trends are unlikely and random unlikely as function of n
trends_threshholds = [];
autocorrelations = [];

sum_of_residual_list = [];

worthy_n = [];
worthy_residuals = [];

for i = 1 : 1 : dimension2
    sum_of_residuals = 0;
    n = all_possible_n(i);
    %n = 11;
    [x_star, r_star] = NOfit(t,Y,n);
    [dim1,dim2] = size(r_star);
    
    mean_of_residuals = mean(abs(r_star));
    standardized_residuals = arrayfun(@(x) (x-mean_of_residuals),abs(r_star));
    
    % It suits me well syntactically, and performance is not in scope of
    % this assignment, nor does it matter here.
    sum_of_residual_list = [sum_of_residual_list norm(r_star, 2)];
    
    m = length(standardized_residuals);
    n_plus = 0;
    n_minus = 0;
    runs = 0;
    current_run = 'none';


    % This calculates the values needed to test for random signs
    for j = 1 : 1 : m
        if strcmp(current_run, 'none')
            if r_star(j) > 0
                n_plus = n_plus + 1;
                runs = runs + 1;
                current_run = 'positive';
                continue;
            end
            if r_star(j) < 0
                n_minus = n_minus + 1;
                runs = runs + 1;
                current_run = 'negative';
                continue;
            end
        end
        
        if strcmp(current_run, 'negative')
            if r_star(j) > 0
                n_plus = n_plus + 1;
                runs = runs + 1;
                current_run = 'positive';
                continue;
            end
            if r_star(j) < 0
                n_minus = n_minus + 1;
                continue;
            end
        end

        if strcmp(current_run, 'positive')
            if r_star(j) > 0
                n_plus = n_plus + 1;
                continue;
            end
            if r_star(j) < 0
                n_minus = n_minus + 1;
                runs = runs + 1;
                current_run = 'negative';
                continue;
            end
        end
    end

    
    % This tests for random signs
    u_mean = (2*n_plus*n_minus / m) + 1;  
    deviation = ((u_mean-1)*(u_mean-2)) / (m-1);
    test_for_random_signs_value = (abs(runs - u_mean)) / (sqrt(deviation));
    random_likely = false;

    if test_for_random_signs_value <= 1.96
        random_likely = true; % With significance level of 0.05 = 5%
    end

    % This tests for correlation
    autocorrelation = 0;
    first_part_of_trend_threshhold = 1 / (sqrt(m-1));
    sum_part_of_threshhold = 0;
    trend_treshhold = 0;
    trend_likely = false; % There should be no trends in the residuals.

    for z = 1 : 1 : m-1
        autocorrelation = autocorrelation + standardized_residuals(z)*standardized_residuals(z+1);
    end

    for h = 1 : 1 : m
        sum_part_of_threshhold = sum_part_of_threshhold + (standardized_residuals(h)^2);
    end

    trend_treshhold = first_part_of_trend_threshhold * sum_part_of_threshhold;

    autocorrelations = [autocorrelations abs(autocorrelation)];
    trends_threshholds = [trends_threshholds trend_treshhold];
    
    if abs(autocorrelation) > trend_treshhold
        trend_likely = true;
    end
    
    if ((~trend_likely) && (~random_likely))
        worthy_n = [worthy_n n];
        worthy_residuals = [worthy_residuals norm(r_star, 2)];
    end
    
    figure(2);
    bar(n, test_for_random_signs_value);
    hold on;
end
figure(1);
hold on;
plot(all_possible_n, sum_of_residual_list);
xlabel('n','Fontsize',14);
ylabel('Residual sizes','Fontsize',14);
title('2 norm of residuals as a function of n');
set(gca,'XTick', 1:2:25);
hold off;

figure(2);
xlim=get(gca,'xlim');
hold on;
plot(xlim,[1.96 1.96], 'r');
xlabel('n','Fontsize',14);
ylabel('Random signs value','Fontsize',14);
title('Significance level = 0.05');
set(gca,'XTick', 1:2:25);
hold off

figure(3);
hold on;
[ax,b,p] = plotyy(all_possible_n,autocorrelations,all_possible_n,trends_threshholds,'bar','plot');
p.LineWidth = 1;
p.Color = 'r';
title('Trends as a function of n');
xlabel('n');
ylabel(ax(1),'Absolute autocorrelation value');
ylabel(ax(2),'Trend threshold');
set(gca,'XTick', 1:2:25);
hold off

figure(4);
hold on;
plot(worthy_n, worthy_residuals);
xlabel('Worthy n','Fontsize',14);
ylabel('Worthy residuals','Fontsize',14);
title('Worthy residuals as function of n');
set(gca,'XTick', 1:2:25);
%set(gca,'YTick', 0:10:150);
%set(gca,'YTick', 0:0.1:20);
hold off;
clearvars;

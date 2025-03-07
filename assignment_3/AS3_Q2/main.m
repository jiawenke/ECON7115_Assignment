
clear
clc
close all

%% parameters
a_x = 0; % lower_bound for x
b_x = 1; % upper_bound for x
a_y = 1; % lower_bound for y
b_y = 3; % upper_bound for y



%% Trapezoid
n_x = 100;
n_y = 100;

function integral_value = trapezoid_2d_integral(a_x, b_x, a_y, b_y, n_x, n_y)
   
    f = @(x, y) x.^2 + y.^2 + 2*x./y;

    h_x = (b_x - a_x) / n_x; 
    h_y = (b_y - a_y) / n_y;

    integral_value = 0;


    for i = 0:n_x
        for j = 0:n_y
            x = a_x + i * h_x;
            y = a_y + j * h_y;
            
            if (i == 0 || i == n_x) && (j == 0 || j == n_y)
                weight = 1;
            elseif (i == 0 || i == n_x) || (j == 0 || j == n_y)
                weight = 2;
            else
                weight = 4;
            end
            
            integral_value = integral_value + weight * f(x, y);
        end
    end

  
    integral_value = integral_value * (h_x * h_y) / 4;
end


trapezoid = trapezoid_2d_integral(a_x, b_x, a_y, b_y, n_x, n_y);
disp(['The trapezoid integral of f(x, y) over [0, 1] x [1, 3] is: ', num2str(trapezoid)]);

%% Simpson's Rule
n_x = 999;
n_y = 999;

function integral_value = simpson_2d_integral(a_x, b_x, a_y, b_y, n_x, n_y)
   
    f = @(x, y) x.^2 + y.^2 + 2*x./y;

    h_x = (b_x - a_x) / (n_x-1); 
    h_y = (b_y - a_y) / (n_y-1);

    integral_value = 0;

    for i = 0:n_x
        for j = 0:n_y
            x = a_x + i * h_x;
            y = a_y + j * h_y;
            
            if mod(i, 2) == 0 && mod(j, 2) == 0
                weight = 4;
            elseif mod(i, 2) == 1 && mod(j, 2) == 1
                weight = 16;
            elseif (mod(i, 2) == 0 && mod(j, 2) == 1) || (mod(i, 2) == 1 && mod(j, 2) == 0)
                weight = 8;
            else
                weight = 0;
            end
            
            integral_value = integral_value + weight * f(x, y);
        end
    end

  
    integral_value = integral_value * (h_x * h_y) / 9;
end



simpson = simpson_2d_integral(a_x, b_x, a_y, b_y, n_x, n_y);
disp(['The Simpson integral of f(x, y) over [0, 1] x [1, 3] is: ', num2str(simpson)]);

%% Monte-Carlo
f = @(x, y) x.^2 + y.^2 + 2*x./y;

% The square of the area
V = (b_x - a_x) * (b_y - a_y);
num_samples = 10000;

% Generate random uniform samples
x_samples = a_x + (b_x - a_x) * rand(num_samples, 1);
y_samples = a_y + (b_y - a_y) * rand(num_samples, 1);

f_values = f(x_samples, y_samples);

integral_value = V * mean(f_values);

disp(['The Monter Carlo integral of f(x, y) over [0, 1] x [1, 3] is: ', num2str(integral_value)]);


%% analytical results to check

f = @(x, y) x.^2 + y.^2 + 2*x./y;

analytical = integral2(f, a_x, b_x, a_y, b_y);

disp(['Analytical integral of f(x, y) over [0, 1] x [1, 3]: ', num2str(analytical)]);

%% References
% Reference: ChatGPT provides the formula of high-dim.
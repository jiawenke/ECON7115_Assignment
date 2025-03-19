clear
close all
clc


%% Parameters

m.beta = 0.95; % discouting rate
m.alpha = 0.33;
m.gamma = 0.5;
m.xi = 0.2;
m.delta = 0.07;


%% Deterministic

m.z = 1;
diff = 1;
tol = 1e-6;
cc = 0;
ww = 0.5;
 a = ones(5, 1);
% initial guess
k_vec = [0.25;0.30;0.35;0.40;0.45];
k_prime_star = ones(5, 1) * 0.1;
%x = zeros(2,1);

x0 = [0.5;0.02];
% iteration
while diff>tol&&cc<10000

for i = 1:5
    ki = k_vec(i);  

    % options = optimoptions('fsolve', 'Display', 'off');
    % [h] = fsolve(@(h) func_l(ki, m, h, a), x0, options);
    % 
    % l_star(i) = h(1);
    % k_prime_star(i) = h(2);

    [l_star(i), k_prime_star(i)] = solve_equations_d(m, ki, k_prime_star(i), a);

    c(i) = (ki^m.alpha * l_star(i)^(1 - m.alpha) - m.delta * ki - k_prime_star(i));
    % c(i) = max((ki^m.alpha * l_star(i)^(1 - m.alpha) - m.delta * ki - k_prime_star(i)),0.01);
    U(i) = (c(i)^(1 - m.gamma)) / (1 - m.gamma) - (l_star(i)^(1 + m.xi)) / (1 + m.xi);

    k_prime = k_prime_star(i);
    V_new(i) = U(i) + func_val(m,a,k_prime);
end

% fit and update a
 A = [ones(5, 1), k_vec, k_vec.^2, k_vec.^3, k_vec.^4];
 a_new = (A' * A) \ (A' * V_new');
 a = ww*a + (1-ww) * a_new;

 diff =  max(abs(a_new - a));
 cc=cc+1;
 
end
disp('Optimal [a_0,...,a_4]:');
disp(a);

%% Stochastic
 
tol = 1e-6;
max_iter = 1000;
a = zeros(5, 1);
k_vec = [0.55; 0.6; 0.65; 0.7; 0.75];
z_vec = [0.8; 1.2];
pi_mat = [
    0.8, 0.2 
    0.4, 0.6];
m.S = length(z_vec);
a = repmat(ones(5, 1), 1, m.S);
k_prime_star_L = ones(5, 1) * 0.1;
k_prime_star_H = ones(5, 1) * 0.1;
diff =1;
cc=1;
ww=0.1;
 
% Define the function to solve the optimal labor and k'
function [l_star, k_prime_star] = solve_equations_s(m, ki, z, k_prime_star_prev, a)
    
    eq1 = @(l) l^(m.xi + m.alpha) - (z * ki^m.alpha * l^(1 - m.alpha) - m.delta * ki - k_prime_star_prev)^(-m.gamma) * (1 - m.alpha) * z * ki^m.alpha;
    options = optimoptions('fsolve', 'Display', 'off');
    l_star = fsolve(eq1, 0.51, options);
 
    
    eq2 = @(k_prime) (z * ki^m.alpha * l_star^(1 - m.alpha) - m.delta * ki - k_prime)^(-m.gamma) - m.beta * (a(2) + 2 * a(3) * k_prime + 3 * a(4) * k_prime^2 + 4 * a(5) * k_prime^3);
    k_prime_star = fsolve(eq2, 0.5, options);
end 

% Define the value function of V(k'*)
function y = func_value(k_prime_star,a)
y = a(1) + a(2)*k_prime_star + a(3)*k_prime_star^2 + a(4)*k_prime_star^3 + a(5)*k_prime_star^4;
end

% iteration
while diff > tol && cc < 1000
    V_new_L = zeros(5, 1);
    V_new_H = zeros(5, 1);
    
    for i = 1:5 
        ki = k_vec(i);
        
        % Solve the optimals
        [l_star_L(i), k_prime_star_L(i)] = solve_equations_s(m, ki, z_vec(1), k_prime_star_L(i), a(:, 1));
        [l_star_H(i), k_prime_star_H(i)] = solve_equations_s(m, ki, z_vec(2), k_prime_star_H(i), a(:, 2));
        
        % Calcualte the consumption
        c_L(i) = ki^m.alpha * l_star_L(i)^(1 - m.alpha) - m.delta * ki - k_prime_star_L(i);
        c_H(i) = ki^m.alpha * l_star_H(i)^(1 - m.alpha) - m.delta * ki - k_prime_star_H(i);
        
        % Calculate the Utility of current period
        U_L(i) = (c_L(i)^(1 - m.gamma)) / (1 - m.gamma) - (l_star_L(i)^(1 + m.xi)) / (1 + m.xi);
        U_H(i) = (c_H(i)^(1 - m.gamma)) / (1 - m.gamma) - (l_star_H(i)^(1 + m.xi)) / (1 + m.xi);
        
        % Calculate the Value of current period
        V_new_L(i) = U_L(i) + m.beta * pi_mat(1, 1)*func_value(k_prime_star_L(i),a(:,1))+ pi_mat(1, 2)*func_value(k_prime_star_L(i),a(:,2));  
        V_new_H(i) = U_H(i) + m.beta * pi_mat(1, 1)*func_value(k_prime_star_H(i),a(:,1))+ pi_mat(1, 2)*func_value(k_prime_star_H(i),a(:,2));
    end 
    
    % Fit and update a
    A = [ones(5, 1), k_vec, k_vec.^2, k_vec.^3, k_vec.^4];
    V_new = [V_new_L, V_new_H];
    a_new = (A' * A) \ (A' * V_new);

    diff = max(abs(a_new(:) - a(:)));
    cc = cc + 1;

    a = ww * a + (1 - ww) * a_new;
end 
 
disp('Optimal a:');
disp(a);

clear
clc
close all

%% Data and Parameters
m.N = 2;
m.R = 4;
m.sigma = 5;
m.theta = m.sigma-1;
m.barA = 1;
m.barB = 1;
m.barL = 1;
m.alpha = 0.1;
m.mu = [0.5;0.5;0.5;0.5];

tau = csvread('data_spatial/tau.csv', 1,1);
m.tau = tau;



%%  eq'm when tariff = 0
tar_1 = 0;
m.tar_2=0;
[w,X,P,L1,welfare] = s_func_eqm_iter(tar_1,m);
disp(welfare)
disp(w)

%% Unilateral Optimal tariff: grid the range
tar_vec = (0:0.02:30)';
T = length(tar_vec);

wel_vec = ones(T,1);
L_mat = ones(T,m.R);

for tt = 1:T
    tar = tar_vec(tt);
    [w,X,P,L1,welfare] = s_func_eqm_iter(tar,m);
    
    wel_vec(tt) = welfare(1);
    L_mat(tt,:) = L1;

    disp(tar)
end

figure(1)
qx = tar_vec;
h = line(qx, wel_vec, 'Color', 'r', 'LineStyle', '-.', 'LineWidth', 2);

% Display results
[max_value, max_index] = max(wel_vec);
tar_1 = qx(max_index);

disp(['Max Welfare ', num2str(max_value)]);
disp(['When Tariff_1 = ', num2str(tar_1), ', Welfare_1 is maximized']);

%% Unilateral Optimal tariff: fminunc
m.tar_2=0;
function neg_welfare = welfare_objective(tar_1, m)
    [~, ~, ~, ~, welfare] = s_func_eqm_iter(tar_1, m);
    neg_welfare = -welfare(1)-welfare(2);
end

initial_tariff1 = 0; 
options = optimoptions('fminunc', 'Display', 'iter', 'TolFun', 1e-6);
[optimal_tariff, min_neg_welfare] = fminunc(@(tar_1) welfare_objective(tar_1, m), initial_tariff1, options);

max_welfare = -min_neg_welfare;
tar_1 = optimal_tariff;
m.tar_1 = optimal_tariff;
disp(['Max Welfare: ', num2str(max_welfare)]);
disp(['When Tariff_1 = ', num2str(optimal_tariff), ', Welfare_1 is maximized']);

%% Nash eq'm

m.tar_1=0;
m.tar_2=0;

options = optimoptions('fminunc', 'Display', 'iter', 'TolFun', 1e-6);

tol = 1e-6;
diff = 1;
cc = 0;

initial_tariff1 = 0;
initial_tariff2 = 0;

while diff > tol && cc < 100
    
    % update country1's tariff
    [optimal_tariff1] = fminunc(@(tar_1) welfare_objective(tar_1, m), initial_tariff1, options);
    m.tar_1 = optimal_tariff1;
    
    % update country2's tariff
    [optimal_tariff2] = fminunc(@(tar_2) welfare_objective2(tar_2, m), initial_tariff2, options);
    m.tar_2 = optimal_tariff2;

    %difference
    diff = max(abs(optimal_tariff1 - initial_tariff1), abs(optimal_tariff2 - initial_tariff2));
    cc = cc + 1;
end

nash_tar = [optimal_tariff1;optimal_tariff2];
disp(nash_tar)

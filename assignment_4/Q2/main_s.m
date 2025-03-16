clear
clc
close all

%% Data and Parameters
m.N = 2;
m.R = 4;
m.sigma = 5;
m.theta = m.sigma-1;
m.barA = ones(m.R,1);
m.barB = ones(m.R,1);
m.barL = ones(m.R,1);
m.alpha = 0.1;
m.mu = [0.5;0.5;0.5;0.5];

tau = csvread('data_spatial/tau.csv', 1,1);
m.tau = tau;

tar_1 = 0;
tar_2 = 0;
tarr_vec = func_tar(m,tar_1,tar_2);

%%  eq'm when tariff = 0


[w,X,P,L,welfare] = func_eqm_iter(tarr_vec,m);
disp(welfare)
disp(w)

%% Unilateral Optimal tariff: grid the range
tar_vec = (0:0.02:0.9)';
T = length(tar_vec);

wel_vec = ones(T,1);
L_mat = ones(T,m.R);

tar_2 = 0;
for tt = 1:T
    tar_1 = tar_vec(tt);

    tarr_vec = func_tar(m,tar_1,tar_2);
    [w,X,P,L,welfare] = func_eqm_iter(tarr_vec,m);
    
    wel_vec(tt) = welfare(1)+welfare(2);
    L_mat(tt,:) = L';

    disp(tar_1)
end

[max_value, max_index] = max(wel_vec);
tar_1_op = tar_vec(max_index);

disp(['Max Welfare ', num2str(max_value)]);
disp(['When Tariff_1 = ', num2str(tar_1_op), ', Welfare_1 is maximized']);

%% Unilateral Optimal tariff: fminunc

function neg_welfare = welfare_objective(m,tar_1,tar_2)
    tarr_vec = func_tar(m,tar_1,tar_2);
    [~, ~, ~, ~, welfare] = func_eqm_iter(tarr_vec, m);
    neg_welfare = -welfare(1)-welfare(2);
end

tar_2=0;
tarr_vec_0 = zeros(4,4);
tarr_vec = func_tar(m,tar_1,tar_2);

options = optimoptions('fminunc', 'Display', 'iter', 'TolFun', 1e-6);
[optimal_tariff, min_neg_welfare] = fminunc(@(tar_1) welfare_objective(m,tar_1,tar_2), 0, options);

max_welfare = -min_neg_welfare;


disp(['Max Welfare: ', num2str(max_welfare)]);
disp(['When Tariff_1 = ', num2str(optimal_tariff), ', Welfare_1 is maximized']);

%% Nash eq'm

tar_1=0;
tar_2=0;

options = optimoptions('fminunc', 'Display', 'iter', 'TolFun', 1e-6);

tol = 1e-6;
diff = 1;
cc = 0;

initial_tariff1 = 0;
initial_tariff2 = 0;

function neg_welfare = welfare_objective_both(m, tar_1, tar_2, country_index)
    tarr_vec = func_tar(m, tar_1, tar_2);
    [~, ~, ~, ~, welfare] = func_eqm_iter(tarr_vec, m);
    
    if country_index == 1
        neg_welfare = -welfare(1) - welfare(2); % country 1
    elseif country_index == 2
        neg_welfare = -welfare(3) - welfare(4); % country 2
    else
        error('Invalid country index. Use 1 or 2.');
    end
end



while diff > tol && cc < 100
    
    % update country1's tariff
    [optimal_tariff1] = fminunc(@(tar_1) welfare_objective_both(m,tar_1,tar_2,1), initial_tariff1, options);
    tar_1 = optimal_tariff1;
    
    % update country2's tariff
    [optimal_tariff2] = fminunc(@(tar_2) welfare_objective_both(m,tar_1,tar_2,2), initial_tariff2, options);
    tar_2 = optimal_tariff2;

    %difference
    diff = max(abs(optimal_tariff1 - initial_tariff1), abs(optimal_tariff2 - initial_tariff2));
    cc = cc + 1;
end

nash_tar = [optimal_tariff1;optimal_tariff2];
disp(nash_tar)

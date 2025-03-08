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

%tar = zeros(m.R,m.R);
%m.tar = tar;




%%  eq'm when tariff = 0
tar = 9;
[w,X,P,L1,welfare] = s_func_eqm_iter(tar,m);
disp(welfare)
disp(w)

%% Unilateral Optimal tariff
tar_vec = (0:50:3800)';
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

%%
tar_2 = 9;
[w,X,P,L1,welfare] = s_func_eqm_iter1(tar_2,m);
disp(welfare)
disp(w)

%% Nash Eqm Tariff

tar_vec = (0:100:5000)';
T = length(tar_vec);

wel_vec = ones(T,1);
L_mat = ones(T,m.R);


for tt = 1:T
    tar_2 = tar_vec(tt);
    [w,X,P,L1,welfare] = s_func_eqm_iter1(tar_2,m);
    
    wel_vec(tt) = welfare(1);
    L_mat(tt,:) = L1;

    disp(tar_2)
end

qx = tar_vec;
h = line(qx, wel_vec, 'Color', 'r', 'LineStyle', '-.', 'LineWidth', 2);

% Display results
[max_value, max_index] = max(wel_vec);
tar_2 = qx(max_index);

disp(['Max Welfare ', num2str(max_value)]);
disp(['When Tariff_2 = ', num2str(tar_2), ', Welfare_2 is maximized']);
